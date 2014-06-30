class Bio::Velvet::Graph::NodedRead
  def adjusted_position(parent_node)
    if @direction == true
      return @offset_from_start_of_node
    elsif @direction == false
      return parent_node.length - @offset_from_start_of_node
    else
      raise "programming error"
    end
  end
end

module Bio
  module AssemblyGraphAlgorithms
    class ContigPrinter
      include Bio::FinishM::Logging

      class AnchoredConnection
        # The identifiers of the probe reads in the velvet assembly graph
        attr_accessor :start_probe_noded_read, :end_probe_noded_read

        # The nodes in the graph that contain the start_probe_read_id and end_probe_read_id
        attr_accessor :start_probe_node, :end_probe_node

        # number of nucleotides between the start of the start probe read and the start of the end of the contig
        attr_accessor :start_probe_contig_offset

        # number of nucleotides until the end of the end probe read in the start of the second contig
        attr_accessor :end_probe_contig_offset

        # Length of the start and end probe sequences
        attr_accessor :start_probe_read_length
        attr_accessor :end_probe_read_length

        # Enumerable of Enumerables of OrientedNode objects, each list of OrientedNode objects
        # corresponds to a path that forms the connection
        attr_accessor :paths
      end

      # Take an Enumerable contigs_and_connections which contains
      # either contigs (as correctly oriented Strings) or an connection (an AnchoredConnection).
      # Cannot handle circular contigs,
      # and it must be in the order contig,connection,contig and start and end with a contig
      def contigs_and_connections_to_string(graph, contigs_and_connections)
        to_return = nil
        last_connection = nil

        contigs_and_connections.each_with_index do |concon, i|
          if last_element.nil?
            raise "Unexpected first element" unless concon.kind_of?(String)
            to_return = concon
          elsif i % 2 == 1
            # Every second element is a connection
            raise "Unexpectedly didn't find AnchoredConnection" unless concon.kind_of?(AnchoredConnection)
            last_connection = concon
          else
            raise "Unexpectedly didn't find String" unless concon.kind_of?(String)
            raise if second_last_element.nil? or last_element.nil?
            to_return = one_connection_between_two_contigs(graph, to_return, last_connection, concon)
          end
        end
        return to_return
      end

      # Given two contigs, return a String representing the new contig. Assumes
      # that there is only 1 path between the two contigs.
      #
      #          ---------->         <--------           start and end probes (ends of probe sequences may not form part of final path). Directions not variable.
      #  --------------------->NNNN------------------->  original sequence to be gapfilled (contig1, NNNN, contig2). Directions not variable
      #      -----------                 ------->        path across the gap. Direction not variable
      #                 \               /
      #                  --------------
      #      ---------->|<-----|----->|--------->        nodes that make up the path (directions and boundaries variable)
      #    stage1|           stage2           |stage3    stages of sequence construction in this method
      def one_connection_between_two_contigs(graph, contig1, anchored_connection, contig2)
        to_return = ''

        log.debug "Working with anchored_connection: #{anchored_connection.inspect}" if log.debug?

        # Stage1 - contig1 before the path begins
        to_return = contig1[0...-(anchored_connection.start_probe_contig_offset)]
        log.debug "After first chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?

        # Stage2 - path sequence, beginning and ending with
        # beginning and ending probes
        begin
          path = nil
          anchored_connection.paths.each do |pat|
            raise "Found multiple paths - can't (yet?) handle this" unless path.nil?
            path = pat
          end

          # Find start index
          begin_node = anchored_connection.start_probe_node
          begin_noded_read = anchored_connection.start_probe_noded_read
          raise if begin_noded_read.nil?
          if begin_noded_read.start_coord != 0
            log.error "Unexpectedly the start of the begin probe not did not form part of the path, possibly indicating misassembly and use of untested code: #{begin_noded_read.inspect}"
            log.error "Anchored connection was: #{anchored_connection.inspect}"
            raise "some kind of error"
            #   to_return += contig1[-(anchored_connection.start_probe_contig_offset)...-(anchored_connection.start_probe_contig_offset+1)]
          end
          offset_of_begin_probe_on_path = begin_noded_read.offset_from_start_of_node

          # Find end index
          end_node = anchored_connection.end_probe_node
          end_noded_read = anchored_connection.end_probe_noded_read
          raise if end_noded_read.nil?
          if end_noded_read.start_coord != 0
            log.error "Unexpectedly the start of the begin probe not did not form part of the path, possibly indicating misassembly and use of untested code: #{end_noded_read.inspect}"
            log.error "Anchored connection was: #{anchored_connection.inspect}"
            raise "some kind of error"
          end
          offset_of_end_node_on_path = end_noded_read.offset_from_start_of_node

          log.debug "Found start index #{offset_of_begin_probe_on_path} and end index #{offset_of_end_node_on_path}" if log.debug?
          path_sequence = path.sequence
          log.debug "Path has a sequence length #{path_sequence.length}" if log.debug?
          log.debug "Returned path sequence will be #{path_sequence.length-offset_of_begin_probe_on_path-offset_of_end_node_on_path} long" if log.debug?
          to_return += path_sequence[offset_of_begin_probe_on_path...-(offset_of_end_node_on_path)]
          log.debug "After path chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?
        end #end stage 2

        # Stage 3
        to_return += contig2[(anchored_connection.end_probe_contig_offset+1)..-1]
        log.debug "After last chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?

        return to_return
      end

      class Variant
        attr_accessor :reference_oriented_node_before_variant, :reference_oriented_node_after_variant, :variation_path

        def to_settable
          [@reference_oriented_node_before_variant.to_settable, @reference_oriented_node_after_variant.to_settable, @variation_path.collect{|onode| onode.to_settable}].flatten
        end
      end

      class Connection
        attr_accessor :reference_path, :variants
      end


      # Given paths between two nodes, return the reference sequence
      # and variants from that reference sequence. paths is an Enumerable
      # where each element is an independent path from the start node to the
      # end node. It is assumed that the start node and the end nodes
      # are the same for each path.
      def two_contigs_and_connection_to_printable_connection(paths)
        # Take the first path as the reference
        reference_path = paths[0]

        #For the moment assume all variants are independent
        all_variants = []
        variant_set = Set.new

        reference_onodes_to_indices = {}
        reference_path.each_with_index do |onode, i|
          reference_onodes_to_indices[onode.to_settable] = i
        end

        # Find variants
        paths.each_with_index do |path, path_i|
          next if path_i==0 #first path is the reference

          # Effectively we are trying to solve an alignment problem here, which is hard.
          # Take the easy route here. This won't well or at all when there is cycles.
          # Assume each common node between the reference and this path 'match up',
          # and the variants are just the bits in between
          current_variant = nil
          previous_reference_onode_index = -1
          path.each_with_index do |onode, onode_i|
            matching_reference_node_index = reference_path[previous_reference_onode_index+1...reference_path.length].find_index(onode)
            if matching_reference_node_index
              matching_reference_node_index += previous_reference_onode_index+1
            end

            log.debug "found matching reference node index #{matching_reference_node_index}" if log.debug?
            log.debug "previous_reference_onode_index #{previous_reference_onode_index}, current_variant: #{current_variant.inspect}" if log.debug?
            if matching_reference_node_index.nil?
              # This node is variant
              if current_variant.nil?
                # new fork. Setup the variant
                current_variant = Variant.new
                raise "not all paths start at the same node!" if previous_reference_onode_index < 0
                current_variant.reference_oriented_node_before_variant = reference_path[previous_reference_onode_index]
                current_variant.variation_path = [onode]
                log.debug "New variant: #{current_variant.inspect}"
              else
                # Building on a current_variant
                current_variant.variation_path << onode
              end
            else
              # Not in a variation (any more?)
              if current_variant.nil?
                # not currently in any variant, and still not.
              else
                # ending a variant
                current_variant.reference_oriented_node_after_variant = onode
                unless variant_set.include?(current_variant.to_settable)
                  all_variants.push current_variant
                  variant_set << current_variant.to_settable
                end
                current_variant = nil
              end
              previous_reference_onode_index = matching_reference_node_index
            end
          end
        end

        to_return = Connection.new
        to_return.reference_path = reference_path
        to_return.variants = all_variants

        return to_return
      end
    end
  end
end

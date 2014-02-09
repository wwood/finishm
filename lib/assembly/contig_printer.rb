

module Bio
  module AssemblyGraphAlgorithms
    class ContigPrinter
      include Bio::FinishM::Logging

      class AnchoredConnection
        # The identifiers of the probe reads in the velvet assembly graph
        attr_accessor :start_probe_read_id, :end_probe_read_id

        # The nodes in the graph that contain the start_probe_read_id and end_probe_read_id
        attr_accessor :start_probe_node, :end_probe_node

        # number of nucleotides between the start of the probe read and the start of the end of the contig
        attr_accessor :start_probe_contig_offset

        # number of nucleotides until the start of the probe read in the start of the second contig
        attr_accessor :end_probe_contig_offset

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
      #TODO: this method will almost certainly fail because it takes no notice of NodedRead.direction
      def one_connection_between_two_contigs(graph,contig1,anchored_connection,contig2)
        to_return = nil

        # add the contig up until the probe read begins
        # then add the bits of the probe that are before the path beginning (if this is non-zero perhaps it indicates a mistake?, or sequencing error in read)
        # add the nucleotides up until the first probe bites in
        # In the velvet graph, this is where the read first starts on the first node
        begin_node = anchored_connection.start_probe_node
        begin_node_read = begin_node.short_reads.find{|noded_read| noded_read.read_id == anchored_connection.start_probe_read_id}
        log.debug "begin noded read: #{begin_node_read.inspect}" if log.debug?
        # class NodedRead
        #   attr_accessor :read_id, :offset_from_start_of_node, :start_coord, :direction
        # end
        first_bite_length = begin_node_read.start_coord - 1
        first_bite_start = contig1.length-anchored_connection.start_probe_contig_offset
        first_bite = contig1[first_bite_start..(first_bite_start+first_bite_length)]
        log.debug "Adding first bite #{first_bite}" if log.debug?
        to_return = contig1[0...first_bite_start] #up until the first bite
        log.debug "Before first bite, sequence is #{to_return}" if log.debug?
        to_return += first_bite #the first bite
        log.debug "Adding first bite sequence `#{first_bite}'" if log.debug?

        # then add the path itself
        path = nil
        anchored_connection.paths.each do |pat|
          raise "Found multiple paths - can't yet handle this" unless path.nil?
          path = pat
        end
        pp path
        raise "path not valid - wrong start node" unless path[0].node == anchored_connection.start_probe_node
        raise "path not valid - wrong end node" unless path[path.length-1].node == anchored_connection.end_probe_node
        whole_seq = path.sequence
        log.debug "Found whole sequence of path #{whole_seq}" if log.debug?
        end_node = anchored_connection.end_probe_node
        end_probe_read = end_node.short_reads.find{|noded_read| noded_read.read_id == anchored_connection.end_probe_read_id}


        # TODO: add some notion of direction to this check
        raise "unhandled 1 node path with confusing start and stop" if end_node == begin_node and
          end_probe_read.offset_from_start_of_node < begin_node_read.offset_from_start_of_node
        path_seq_start = begin_node_read.offset_from_start_of_node
        path_seq_end = whole_seq.length - end_probe_read.offset_from_start_of_node
        log.debug "path_seq_start=#{path_seq_start}, path_seq_end=#{path_seq_end}, end_probe_read.offset_from_start_of_node=#{end_probe_read.offset_from_start_of_node}"
        to_return += whole_seq[path_seq_start...path_seq_end]
        log.debug "after adding the path's sequence, now have #{to_return}" if log.debug?

        # then add the bits after
        second_bite_start = anchored_connection.end_probe_contig_offset-end_probe_read.offset_from_start_of_node
        to_return += contig2[second_bite_start...anchored_connection.end_probe_contig_offset]
        log.debug "after adding second bite, sequence is #{to_return}" if log.debug?

        # then add the second contig, which will be a chopped in at both ends
        to_return += contig2[anchored_connection.end_probe_contig_offset...contig2.length]
        log.debug "after adding second contig, sequence is #{to_return}" if log.debug?

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

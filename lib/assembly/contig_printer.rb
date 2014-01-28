

module Bio
  module AssemblyGraphAlgorithms
    class ContigPrinter
      include Bio::FinishM::Logging


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

      class Variant
        attr_accessor :reference_oriented_node_before_variant, :reference_oriented_node_after_variant, :variation_path

        def to_settable
          [@reference_oriented_node_before_variant.to_settable, @reference_oriented_node_after_variant.to_settable, @variation_path.collect{|onode| onode.to_settable}].flatten
        end
      end

      class Connection
        attr_accessor :reference_path, :variants
      end
    end
  end
end

require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms

    # Represents a set of trails, and whether or not circularity has been detected.
    class TrailSet
      attr_accessor :trails
      attr_accessor :circular_paths_detected
      include Enumerable

      def each
        @trails.each{|t| yield t}
      end
    end

    class AcyclicConnectionFinder
      include Bio::FinishM::Logging

      # Find trails between two oriented nodes, both facing the same way along the path.
      #
      # Options:
      # * :recoherence_kmer: use a longer kmer to help de-bubble and de-cicularise (default don't use this)
      # * :sequences: Bio::Velvet::Sequence object holding sequences of nodes within leash length
      def find_trails_between_nodes(graph, initial_oriented_node, terminal_oriented_node, leash_length, options={})

        #TODO: this is now implemented in the finishm_graph object - just get it from there
        initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        initial_path.add_oriented_node initial_oriented_node

        finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
        return finder.find_all_connections_between_two_nodes(
          graph, initial_path, terminal_oriented_node, leash_length, options[:recoherence_kmer], options[:sequences], options
          )
      end

      # Algorithms like SingleCoherentWanderer#wander give an overly short
      # base pair distance between two probes, because the length of the node
      # containing the probe at either end is not included in the calculation.
      #
      # Return the calibrated distance i.e. the true base pair distance between
      # the start of each node pair. Returned is the given distance plus the
      # distance between the start of each probe and the end of the containing
      # node.
      def calibrate_distance_accounting_for_probes(finishm_graph, probe1_index, probe2_index, distance)
        read1 = finishm_graph.probe_node_reads[probe1_index]
        read2 = finishm_graph.probe_node_reads[probe2_index]
        probe_node1 = finishm_graph.probe_nodes[probe1_index]
        probe_node2 = finishm_graph.probe_nodes[probe2_index]

        # If the start and end nodes are the same, that's a special case:
        if finishm_graph.probe_nodes[probe1_index].node_id == finishm_graph.probe_nodes[probe2_index].node_id
          if (read1.direction == true and read2.direction == false) or
            (read1.direction == false and read2.direction == true)
            return probe_node1.length - read1.offset_from_start_of_node - read2.offset_from_start_of_node - finishm_graph.graph.hash_length
          else
            raise "Programming error: to connect within a single contig two probes must have opposite directions: found #{read1.direction} and #{read2.direction}"
          end
        else
          # Usual case - start and end nodes are different nodes
          to_return = distance
          # add the first probe side
          to_return += probe_node1.length-read1.offset_from_start_of_node-finishm_graph.graph.hash_length
          # add the second probe side
          to_return += probe_node2.length-read1.offset_from_start_of_node-finishm_graph.graph.hash_length
          return to_return
        end
      end
    end
  end
end

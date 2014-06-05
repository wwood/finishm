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

        if options[:recoherence_kmer]
          finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
          return finder.find_all_connections_between_two_nodes(
            graph, initial_path, terminal_oriented_node, leash_length, options[:recoherence_kmer], options[:sequences]
            )
        else
          return Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new.find_all_connections_between_two_nodes(
            graph, initial_path, terminal_oriented_node, leash_length
            )
        end
      end

      # Takes a set of probes, and try to see which ones connect together within the leash length.
      #
      # However this method does not determine the paths between each pair of nodes,
      # it merely shows that a path exists.
      #
      # Return a hash repesenting the connections found:
      #
      #     [probe_index1, probe_index2] => min_distance
      #
      # probe_index1 will always be less than probe_index2.
      #
      # TODO: the method used here is probably sub-optimal. It might be better implemented
      # with Dijkstra's shortest path algorithm, or perhaps Johnson's algorithm if that can
      # be modified to only find distances between a subset of nodes in the graph.
      def depth_first_search_with_leash(finishm_graph, leash_length)
        to_return = {}

        # Take the probes and make them all into finishing nodes
        finishing_nodes = []
        finishm_graph.probe_nodes.each_with_index do |probe_node, probe_node_index|
          direction = finishm_graph.probe_node_directions[probe_node_index]
          if direction == true
            finishing_nodes.push [probe_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST]
          else
            finishing_nodes.push [probe_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]
          end
        end

        # Do a depth first search from each probed node in the graph
        finishm_graph.probe_nodes.each_with_index do |probe_node, probe_node_index|
          # Do a depth first search starting from this node. Go all the way to the leash length,
          # and then search to see if any of the other nodes have been come across
          log.info "Exploring from probe node \##{probe_node_index+1} (node #{probe_node.node_id}/#{finishm_graph.probe_node_directions[probe_node_index] })"
          known_nodes_set = Set.new #set of oriented nodes (but not distances - these are stored separately)
          depth_first_search_stack = DS::Stack.new
          initial = finishm_graph.initial_path_from_probe(probe_node_index)
          if initial.nil?
            log.warn "Unable to start searching from probe \##{probe_node_index+1}, because it was no found in the graph. Skipping."
            next
          end
          initial_distanced = DistancedOrientedNode.new
          initial_distanced.node = initial.last.node
          initial_distanced.first_side = initial.last.first_side
          initial_distanced.distance = 0

          minimum_node_distances = {}

          depth_first_search_stack.push initial_distanced
          while distanced_node = depth_first_search_stack.pop
            if known_nodes_set.include?(distanced_node.to_settable) and
              distanced_node.distance >= minimum_node_distances[distanced_node.to_settable]
              # This node has already been explored, and no shorter path has been found. Go no further.
              next
            else
              known_nodes_set << distanced_node.to_settable
            end
            minimum_node_distances[distanced_node.to_settable] = distanced_node.distance

            if distanced_node.distance <= leash_length
              # Still within the leash. Push into the stack all the current node's neighbours in the graph
              trail = Bio::Velvet::Graph::OrientedNodeTrail.new
              trail.add_node distanced_node.node, distanced_node.first_side
              trail.neighbours_of_last_node(finishm_graph.graph).each do |oriented_neighbour|
                d = DistancedOrientedNode.new
                d.node = oriented_neighbour.node
                d.first_side = oriented_neighbour.first_side
                d.distance = distanced_node.distance + distanced_node.node.length_alone
                depth_first_search_stack.push d
              end
            else
              # we are beyond the leash, go no further
            end
          end

          # Finished depth first search. Now we have a set of nodes and distances.
          # Check if the nodes
          finishm_graph.probe_nodes.each_with_index do |node, i|
            next if i <= probe_node_index # only return the 'upper triangle' of the distance matrices

            finish = finishing_nodes[i]

            # Disallow hitting probes if they are both on the same node, but
            # not facing each other.
            # -->    <--- OK
            # <--    --> not ok
            # <--    <-- not ok
            # -->    --> not ok
            #binding.pry
            #TODO: deal with the bug where 2 probes are on the same node

            if minimum_node_distances.key?(finish)
              min_distance = minimum_node_distances[finish]
              log.info "Found a connection between probes #{probe_node_index} and #{i}, distance: #{min_distance}"
              to_return[[probe_node_index, i]] = min_distance
            end
          end
        end
        return to_return
      end

      # An oriented node some distance from the origin of exploration
      class DistancedOrientedNode
        attr_accessor :node, :first_side, :distance

        # Using Set object, often we want two separate objects to be considered equal even if
        # they are distinct objects
        def to_settable
          [@node.node_id, @first_side]
        end

        # Which side of the node is not first?
        def second_side
          @first_side == OrientedNodeTrail::START_IS_FIRST ?
            OrientedNodeTrail::END_IS_FIRST :
            OrientedNodeTrail::START_IS_FIRST
        end
      end
    end
  end
end

require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class GraphExplorer
      # Return all paths that emenate from a given node, in the graph
      def explore_from_node(graph, initial_path, leash_length)
        # Do a simple depth first search, forking at each node. Vanilla graph traversal.
        depth_first_search_stack = DS::Stack.new
        first_path = ExplorationPath.new initial_path
        depth_first_search_stack.push first_path
        found_paths = []
        # While there's more paths to explore
        while current_path = depth_first_search_stack.pop
          last = current_path.path.last
          if !leash_length.nil? and current_path.path.length_in_bp > leash_length
            current_path.termination_type = 'Leashed'
            found_paths.push current_path
          else
            neighbours = current_path.path.neighbours_of_last_node(graph)
            if neighbours.empty?
              current_path.termination_type = 'Dead end / coverage'
              found_paths.push current_path
            else
              neighbours_to_add = []
              neighbours.each do |oriented_neighbour|
                # Test for loops, I'm only interested in acyclic paths for the moment
                if current_path.include?(oriented_neighbour)
                  #loop found, terminate path
                  new_path = current_path.copy
                  new_path.add_node oriented_neighbour
                  new_path.termination_type = 'Loop'
                  found_paths.push new_path
                else
                  neighbours_to_add.push oriented_neighbour
                end
              end
              neighbours_to_add.each_with_index do |oriented_neighbour, i|
                # If the last neighbour is being added here, reuse the path
                next_path = nil
                if i == neighbours_to_add.length-1
                  next_path = current_path
                else
                  next_path = current_path.copy
                end
                next_path.add_node oriented_neighbour
                depth_first_search_stack.push next_path
              end
            end
          end
        end

        return found_paths
      end

      class ExplorationPath
        attr_accessor :path, :set_of_nodes, :termination_type

        def initialize(path)
          @path = path
          @set_of_nodes = Set.new path.collect{|n| n.to_settable}
        end

        def include?(oriented_node)
          @set_of_nodes.include?(oriented_node.to_settable)
        end

        def add_node(onode)
          path.add_oriented_node onode
          @set_of_nodes << onode.to_settable
        end

        def copy
          anew = ExplorationPath.new @path.copy
          return anew
        end

        def to_s
          @path.collect{|on| on.node_id}.join(',')
        end
      end
    end
  end
end

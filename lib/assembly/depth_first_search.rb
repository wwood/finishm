require 'ds'
require 'set'

module Bio
  module Velvet
    class Graph
      # Do a depth first search starting at this oriented node,
      # yielding an OrientedNodeTrail at each new node encountered.
      # The new node is the last node in the yielded trail. The result
      # of the yield tells whether this method whether to abandon
      # searching further from this point (false)
      # or keep going (true).
      def depth_first_search(oriented_node)
        log = Bio::Log::LoggerPlus['finishm']
        discovered_oriented_nodes = Set.new

        # Traverse through the graph, yielding at each new node
        current_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        current_path.add_oriented_node oriented_node

        stack = DS::Stack.new
        stack.push current_path

        # While there is more on the stack
        while current_path = stack.pop
          log.debug "Perhaps #{current_path.last}?" if log.debug?
          if discovered_oriented_nodes.include?(path_to_searchable(current_path))
            # Already seen this node, do nothing with it
            log.debug "Skipping #{current_path.last} since that has already been seen" if log.debug?
            next
          else
            log.debug "That's a new node, #{current_path.last}" if log.debug?
            # Found a new node for the user to play with
            discovered_oriented_nodes << path_to_searchable(current_path)

            continue = yield current_path

            # prep for next time if required.
            if continue
              # Sort node IDs to simplify testing
              next_nodes = current_path.neighbours_of_last_node(self).sort{|n1, n2|
                -(n1.node.node_id <=> n2.node.node_id)
              }
              next_nodes.each do |n|
                path = current_path.copy
                path.add_oriented_node n
                stack.push path
              end
            end
          end
        end
      end

      private
      # Set#include? doesn't pick up when the same OrientedNode is picked
      # up twice independently, I don't think. So convert to an array first
      def path_to_searchable(path)
        last = path.last
        return [last.node.node_id, last.first_side]
      end
    end
  end
end

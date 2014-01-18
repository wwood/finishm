require 'set'
require 'ds'

module DS
  # Like DS::PriorityQueue except give the ability to define how priority is given
  class AnyPriorityQueue < PriorityQueue
    #Create new priority queue. Internaly uses heap to store elements.
    def initialize
      @store = BinaryHeap.new {|parent,child| yield parent.key, child.key}
    end
  end
end

module Bio
  module AssemblyGraphAlgorithms
    class BubblyAssembler
      # Starting at a node within a graph, walk through the graph
      # accepting forks, so long as the fork paths converge within some finite
      # length in the graph (the leash length, measured in number of base pairs).
      #
      # Return an Array of Path arrays, a MetaPath, where each path array are the different paths
      # that can be taken at each fork point
      def assemble_from_node(velvet_graph, leash_length, starting_path)
        currently_in_bubble = false

        metapath = MetaPath.new
        metapath << [starting_path.copy]

        # Priority queue to determine which path in the graph to explore next
        # Prioritise things that have lesser numbers, not greater numbers as is default
        queue = DS::AnyPriorityQueue.new {|a,b| a<=b}

        while true
          if currently_in_bubble == false
            while oriented_neighbours = metapath.last[0].neighbours_of_last_node(velvet_graph)
              if oriented_neighbours.empty?
                # This is just a straight out dead end, and we can go no further.
                metapath << start_path
                return metapath
              elsif oriented_neighbours.length == 1
                # Linear thing here, just keep moving forward
                start_path.add_oriented_node oriented_neighbours[0]
              else
                # Reached a fork in the graph here, the point of this algorithm, really.
                currently_in_bubble = true
                last_bubble_node = metapath.last[0].last #last node in the path that was
                bubble_paths = oriented_neighbours.collect do |oneigh|
                  path = metapath.last[0].copy
                  path.add_oriented_node oneigh
                  path
                end
                # queue up all the possibilities
                bubble_paths.each{|bp| queue.enqueue bp, 0}
                break #break out of linear path mode while loop
              end
            end

          else
            # We are in a bubble. Go get some.
            raise DOSTUFF
          end
        end
      end

      class MetaPath < Array
        def last
          self[self.length-1]
        end
      end
    end
  end
end

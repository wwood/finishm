require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class PathsBetweenNodesFinder
      include Bio::FinishM::Logging


      def find_all_connections_between_two_nodes(graph, initial_path, terminal_node, leash_length)
        problems = find_all_problems(graph, initial_path, terminal_node, leash_length)
        return find_paths_from_problems(problems)
      end

      def find_all_problems(graph, initial_path, terminal_node, leash_length)
        # setup dynamic programming cache
        problems = {}

        # setup stack to keep track of initial nodes
        stack = DS::Stack.new
        stack.push initial_path.copy

        push_next_neighbours = lambda do |current_path, last|
          next_nodes = current_path.neighbours_of_last_node(graph)

          #TODO: not neccessary to copy all paths, can just continue one of them

          next_nodes.each do |n|
            path = current_path.copy
            path.add_oriented_node n
            stack.push path
          end
        end

        while current_path = stack.pop
          # Have we solved this before? If so, add this path to that solved problem.
          last = current_path.last
          path_length = current_path.length_in_bp

          if problems[last]
            prob = problems[last]
            prob.known_paths.push current_path

            # If a lesser min distance is found, then we need to start exploring from the
            # current place again
            if path_length < prob.min_distance
              prob.min_distance = path_length
              push_next_neighbours.call current_path, last
            end
          elsif !leash_length.nil? and path_length > leash_length
            # we are past the leash length, give up
          else
            # new problem being solved here
            problem = DynamicProgrammingProblem.new
            problem.min_distance = path_length
            problem.known_paths.push current_path.copy
            problems[last] = problem

            # explore the forward neighbours
            push_next_neighbours.call current_path, last
          end
        end

        return problems
      end

      def find_paths_from_problems(problems)
        stack = DS::Stack.new

        # There should be at least one problem that finishes at the terminal node
        raise "not implemented"
      end


      class DynamicProgrammingProblem
        attr_accessor :min_distance, :known_paths

        def initialize
          @known_paths = []
        end
      end
    end
  end
end
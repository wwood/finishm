require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class PathsBetweenNodesFinder
      include Bio::FinishM::Logging


      def find_all_connections_between_two_nodes(graph, initial_path, terminal_node, leash_length)
        problems = find_all_problems(graph, initial_path, terminal_node, leash_length)
        #pp problems
        problems.each do |key, dynprob|
          dynprob.remove_duplication_in_known_paths!
        end
        puts problems.collect{|key, dynprob|
          [
            key[0],
            dynprob.min_distance,
            dynprob.known_paths.collect{|path| path.to_short_s}.join(' = ')
          ].join(' ')
        }.join("\n")

        #exit
        return find_paths_from_problems(problems, terminal_node)
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
          log.debug "considering #{current_path}"
          log.debug "considering last: #{last}"

          if last == terminal_node
          log.debug "last is terminal"
            # we are done, go no further
            # the problem of getting to the last node is already recorded, so don't need to do anything
            problems[last.to_settable] ||= DynamicProgrammingProblem.new
            problems[last.to_settable].known_paths ||= []
            problems[last.to_settable].known_paths.push current_path

          elsif problems[last.to_settable]
            prob = problems[last.to_settable]
            prob.known_paths.push current_path

            # If a lesser min distance is found, then we need to start exploring from the
            # current place again
            if path_length < prob.min_distance
              log.debug "Found a node with min_distance greater than path length.."
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
            problems[last.to_settable] = problem

            # explore the forward neighbours
            push_next_neighbours.call current_path, last
          end
        end
        puts problems.collect{|key, dynprob|
          [
            key[0],
            dynprob.min_distance,
            dynprob.known_paths.collect{|path| path.to_short_s}.join(' = ')
          ].join(' ')
        }.join("\n")

        return problems
      end

      def find_paths_from_problems(problems, terminal_node)
        stack = DS::Stack.new

        # push all solutions to the "ending in the final node" solutions to the stack
        overall_solution = problems[terminal_node.to_settable]
        pp overall_solution
        pp terminal_node
        return [] if overall_solution.nil? # if there is no solutions to the overall problem then there is no solution at all
        stack.push [
          overall_solution.known_paths[0].to_a,
          []
        ]
        all_paths = []
        pp stack

        while path_halves = stack.pop
          log.debug path_halves.collect{|half| half.collect{|onode| onode.node.node_id}.join(',')}.join(' and ')
          first_half = path_halves[0]
          second_half = path_halves[1]
          if first_half.length == 0
            # If we've tracked all the way to the beginning
            all_paths.push second_half
          else
            last = first_half.last
            if second_half.include?(last)
              # Ignore - this is a cycle, which rarely happens
            else
              paths_to_last = problems[last.to_settable].known_paths
              paths_to_last.each do |path|
                to_push = [path[0...(path.length-1)],[last,second_half].flatten]
                log.debug "Pushing #{to_push.collect{|half| half.collect{|onode| onode.node.node_id}.join(',')}.join(' and ')}"
                stack.push to_push
              end
            end
          end
        end
        pp all_paths

        return all_paths
      end


      class DynamicProgrammingProblem
        attr_accessor :min_distance, :known_paths

        def initialize
          @known_paths = []
        end

        # With leash length considerations, sometimes we can get multiple paths leading to
        # duplications ie 1 -> 2 -> 3 -> 4, also 1->5->3->4 - if length of node 5 is less than length
        # of node 2, then there'll be 2 paths attached
        def remove_duplication_in_known_paths!
          second_to_last_node_observations = Set.new
          @known_paths.select! do |path|
            if path.length > 1
              second_to_last = path[path.length-2]
              if second_to_last_node_observations.include?(second_to_last.to_settable)
                false
              else
                second_to_last_node_observations << second_to_last.to_settable
                true
              end
            else
              true #keep all short paths
            end
          end
        end
      end
    end
  end
end

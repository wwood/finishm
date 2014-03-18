require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class SingleCoherentPathsBetweenNodesFinder
      include Bio::FinishM::Logging


      # Find all paths between the initial and terminal node in the graph.
      # Don't search in the graph when the distance in base pairs exceeds the leash length.
      # Recohere reads (singled ended only) in an attempt to remove bubbles.
      def find_all_connections_between_two_nodes(graph, initial_path, terminal_oriented_node, leash_length, recoherence_kmer)
        problems = find_all_problems(graph, initial_path, terminal_oriented_node, leash_length, recoherence_kmer)
        #pp problems
        problems.each do |key, dynprob|
          dynprob.remove_duplication_in_known_paths!
        end
#         puts problems.collect{|key, dynprob|
#           [
#             key[0],
#             dynprob.min_distance,
#             dynprob.known_paths.collect{|path| path.to_short_s}.join(' = ')
#           ].join(' ')
#         }.join("\n")

        #exit
        paths = find_paths_from_problems(problems, terminal_oriented_node)
        return paths
      end

      def find_all_problems(graph, initial_path, terminal_node, leash_length)
        # setup dynamic programming cache
        problems = {}

        # setup stack to keep track of initial nodes
        stack = DS::Stack.new
        stack.push initial_path.copy

        push_next_neighbours = lambda do |current_path|
          next_nodes = current_path.neighbours_of_last_node(graph)
          log.debug "Pushing #{next_nodes.length} new neighbours of #{current_path.last}" if log.debug?
          #TODO: not neccessary to copy all paths, can just continue one of them
          next_nodes.each do |n|
            log.debug "Pushing neighbour to stack: #{n}" if log.debug?
            path = current_path.copy
            path.add_oriented_node n
            stack.push path
          end
        end

        path_to_settable = lambda do |path|
          cumulative_length = 0
          trail = path.trail
          i = trail.length - 1
          while i >= 0 and cumulative_length < max_bases_between_read_starts
            cumulative_length += trail[i].length_alone
            i -= 1
          end
          # 'Return' an array made up of the settables
          trail[i..-1].collect{|t| t.to_settable}.flatten
        end

        while current_path = stack.pop
          # Have we solved this before? If so, add this path to that solved problem.
          set_key = path_to_settable.call current_path
          path_length = current_path.length_in_bp
          log.debug "considering #{current_path}, path length: #{path_length}" if log.debug?
          log.debug "considering last: #{last}" if log.debug?

          # Unless the path validates, forget it.
          next if !validate_last_node_of_path_by_recoherence(current_path, recoherence_kmer)

          if current_path.last == terminal_node
            log.debug "last is terminal" if log.debug?
            problems[set_key] ||= DynamicProgrammingProblem.new
            problems[set_key].known_paths ||= []
            problems[set_key].known_paths.push current_path

          elsif problems[set_key]
            log.debug "Already seen this problem" if log.debug?
            prob = problems[set_key]
            prob.known_paths.push current_path

            # If a lesser min distance is found, then we need to start exploring from the
            # current place again
            if path_length < prob.min_distance
              log.debug "Found a node with min_distance greater than path length.." if log.debug?
              prob.min_distance = path_length
              push_next_neighbours.call current_path
            end
          elsif !leash_length.nil? and path_length > leash_length
            # we are past the leash length, give up
            log.debug "Past leash length, giving up" if log.debug?
          else
            log.debug "New dynamic problem being solved" if log.debug?
            # new problem being solved here
            problem = DynamicProgrammingProblem.new
            problem.min_distance = path_length
            problem.known_paths.push current_path.copy
            problems[set_key] = problem

            # explore the forward neighbours
            push_next_neighbours.call current_path
          end
          log.debug "Stack size: #{stack.size}" if log.debug?
        end

#         puts problems.collect{|key, dynprob|
#           [
#             key[0],
#             key[1],
#             dynprob.min_distance,
#             dynprob.known_paths.collect{|path| path.to_short_s}.join(' = ')
#           ].join(' ')
#         }.join("\n")

        return problems
      end

      # Given an OrientedNodeTrail, and an expected number of
      def validate_last_node_of_path_by_recoherence(path, recoherence_kmer)
        #not possible to fail on a 1 or 2 node path, by debruijn graph definition.
        return true if path.length < 3

        # Walk backwards along the path from the 2nd last node,
        # collecting nodes until the length in bp of the nodes is > recoherence_kmer
        collected_nodes = []
        length_of_nodes = lambda do |nodes|
          nodes.reduce(0) do |sum, node|
            sum += node.node.length_alone
          end
        end
        i = path.trail.length-2
        while length_of_nodes.call(collected_nodes) < recoherence_kmer and i >= 0
          collected_nodes.push path.trail[i]
          i -= 1
        end
        log.debug "validate: Collected nodes: #{collected_nodes}" if log.debug?

        # There should be at least 1 read that spans the collected nodes and the last node
        # The trail validates iff the above statement is true.
        #TODO: there's a possible 'bug' here in that there's garauntee that the read is overlays the
        # nodes in a consecutive and gapless manner. But I suspect that is unlikely to be a problem in practice.
        possible_reads = path.trail[-1].node.short_reads.collect{|nr| nr.read_id}
        log.debug "validate: Initial short reads: #{possible_reads.join(',')}" if log.debug?
        collected_nodes.each do |node|
          current_set = Set.new node.node.short_reads.collect{|nr| nr.read_id}
          possible_reads.select! do |r|
            current_set.include? r
          end
          return false if possible_reads.empty?
        end
        return true
      end

      def find_paths_from_problems(problems, terminal_node)
        stack = DS::Stack.new

        to_return = Bio::AssemblyGraphAlgorithms::TrailSet.new
        to_return.circular_paths_detected = false
        to_return.trails = []

        # push all solutions to the "ending in the final node" solutions to the stack
        overall_solution = problems[terminal_node.to_settable]
        return to_return if overall_solution.nil? # if there is no solutions to the overall problem then there is no solution at all
        stack.push [
          overall_solution.known_paths[0].to_a,
          []
        ]
        all_paths = []

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
              #TODO: circular paths should be dealt with in some manner. Really no simple solution, however,
              # particularly when there is more than one connected circuit detected.
              log.warn "Linking path(s) detected, but cycle also detected. Giving up on this link."
              to_return.circular_paths_detected = true
              return to_return
            else
              paths_to_last = problems[last.to_settable].known_paths
              paths_to_last.each do |path|
                to_push = [path[0...(path.length-1)],[last,second_half].flatten]
                log.debug "Pushing #{to_push.collect{|half| half.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') }" if log.debug?
                stack.push to_push
              end
            end
          end
        end

        to_return.trails = all_paths
        return to_return
      end


      class DynamicProgrammingProblem
        attr_accessor :min_distance, :known_paths

        def initialize
          @known_paths = []
        end

        # With leash length considerations, sometimes we can get multiple paths leading to
        # duplications ie 1 -> 2 -> 3 -> 4, also 1->5->3->4 - if length of node 5 is less than length
        # of node 2, then there'll be 2 paths attached to 4.
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

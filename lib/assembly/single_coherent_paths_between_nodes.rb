require 'ds'
require 'set'

class DS::Queue
  alias_method :size, :length
end

module Bio
  module AssemblyGraphAlgorithms
    class SingleCoherentPathsBetweenNodesFinder
      include Bio::FinishM::Logging

      SINGLE_BASE_REVCOM = {
        'A'=>'T',
        'T'=>'A',
        'G'=>'C',
        'C'=>'G',
        }

      # Find all paths between the initial and terminal node in the graph.
      # Don't search in the graph when the distance in base pairs exceeds the leash length.
      # Recohere reads (singled ended only) in an attempt to remove bubbles.
      def find_all_connections_between_two_nodes(graph, initial_path, terminal_oriented_node,
        leash_length, recoherence_kmer, sequence_hash)

        problems = find_all_problems(graph, initial_path, terminal_oriented_node, leash_length, recoherence_kmer, sequence_hash)
        problems.each do |key, dynprob|
          dynprob.remove_duplication_in_known_paths!
        end

        paths = find_paths_from_problems(problems, recoherence_kmer)
        return paths
      end

      def find_all_problems(graph, initial_path, terminal_node, leash_length, recoherence_kmer, sequence_hash)
        # setup dynamic programming cache
        problems = ProblemSet.new

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

        current_oriented_node_trail = Bio::Velvet::Graph::OrientedNodeTrail.new

        while current_path = stack.pop
          path_length = current_path.length_in_bp
          log.debug "considering #{current_path}, path length: #{path_length}" if log.debug?

          # Have we solved this before? If so, add this path to that solved problem.
          set_key = path_to_settable current_path, recoherence_kmer
          log.debug "Set key is #{set_key}"

          # Unless the path validates, forget it.
          if !validate_last_node_of_path_by_recoherence(current_path, recoherence_kmer, sequence_hash)
            log.debug "Path did not validate, skipping" if log.debug?
            next
          elsif log.debug?
            log.debug "Path validates"
          end

          if current_path.last == terminal_node
            log.debug "last is terminal" if log.debug?
            problems[set_key] ||= DynamicProgrammingProblem.new
            problems[set_key].known_paths ||= []
            problems[set_key].known_paths.push current_path

            problems.terminal_node_keys ||= Set.new
            problems.terminal_node_keys << set_key

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

        return problems
      end

      def path_to_settable(path, recoherence_kmer)
        log.debug "Making settable a path: #{path}" if log.debug?
        return array_trail_to_settable(path.trail, recoherence_kmer)
      end

      def array_trail_to_settable(trail, recoherence_kmer)
        cumulative_length = 0
        i = trail.length - 1
        while i >= 0 and cumulative_length < recoherence_kmer
          cumulative_length += trail[i].node.length_alone
          i -= 1
        end
        i += 1
        # 'Return' an array made up of the settables
        to_return = trail[i..-1].collect{|t| t.to_settable}.flatten
        log.debug "'Returning' settable version of path: #{to_return}" if log.debug?
        to_return
      end

      # Given an OrientedNodeTrail, and an expected number of
      def validate_last_node_of_path_by_recoherence(path, recoherence_kmer, sequence_hash)
        #not possible to fail on a 1 or 2 node path, by debruijn graph definition.
        return true if path.length < 3

        # Walk backwards along the path from the 2nd last node,
        # collecting nodes until the length in bp of the nodes is > recoherence_kmer
        collected_nodes = []
        length_of_nodes = lambda do |nodes|
          if nodes.empty?
            0
          else
            hash_offset = nodes[0].node.parent_graph.hash_length-1
            nodes.reduce(hash_offset) do |sum, node|
              sum += node.node.length_alone
            end
          end
        end
        i = path.length-2
        while i >= 0
          collected_nodes.push path.trail[i]
          i -= 1
          break if length_of_nodes.call(collected_nodes) >= recoherence_kmer
        end
        log.debug "validate: Collected nodes: #{collected_nodes}" if log.debug?

        # There should be at least 1 read that spans the collected nodes and the last node
        # The trail validates iff the above statement is true.
        #TODO: there's a possible 'bug' here in that there's garauntee that the read is overlays the
        # nodes in a consecutive and gapless manner. But I suspect that is unlikely to be a problem in practice.
        final_node = path.trail[-1].node
        possible_reads = final_node.short_reads.collect{|nr| nr.read_id}
        log.debug "validate starting from #{final_node.node_id}: Initial short reads: #{possible_reads.join(',') }" if log.debug?
        collected_nodes.each do |node|
          log.debug "Validating node #{node}"
          current_set = Set.new node.node.short_reads.collect{|nr| nr.read_id}
          possible_reads.select! do |r|
            current_set.include? r
          end
          if possible_reads.empty?
            log.debug "First line validation failed, now detecting sub-kmer sequence overlap" if log.debug?
            trail_to_validate = path.trail[i+1..-1]
            return sub_kmer_sequence_overlap?(trail_to_validate, sequence_hash)
          end
        end
        return true
      end

      # Is there overlap across the given nodes, even if the overlap
      # does not include an entire kmer?
      # nodes: an OrientedNodeTrail. To validate, there must be at least 1 read that spans all of these nodes
      # sequence_hash: Bio::Velvet::Sequence object with the sequences from the reads in the nodes
      def sub_kmer_sequence_overlap?(nodes, sequence_hash)
        raise if nodes.length < 3 #should not get here - this is taken care of above
        log.debug "validating by sub-kmer sequence overlap: #{nodes}" if log.debug?

        # Only reads that are in the second last node are possible, by de-bruijn graph definition.
        candidate_noded_reads = nodes[-2].node.short_reads
        middle_nodes_length = nodes[1..-2].reduce(0){|sum, n| sum += n.node.length}+
          +nodes[0].node.parent_graph.hash_length-1
        log.debug "Found middle nodes length #{middle_nodes_length}" if log.debug?

        candidate_noded_reads.each do |read|
          # Ignore reads that don't come in at the start of the node
          log.debug "Considering read #{read.inspect}" if log.debug?
          if read.offset_from_start_of_node != 0
            log.debug "Read doesn't start at beginning of node, skipping" if log.debug?
            next
          else
            seq = sequence_hash[read.read_id]
            raise "No sequence stored for #{read.read_id}, programming fail." if seq.nil?

            if read.start_coord == 0
              log.debug "start_coord Insufficient length of read" if log.debug?
              next
            elsif seq.length-read.start_coord-middle_nodes_length < 1
              log.debug "other_side Insufficient length of read" if log.debug?
              next
            end
            binding.pry

            # Now ensure that the sequence matches correctly
            # left base, the base from the first node
            first_node = nodes[0].node
            left_base = !(read.direction ^ nodes[-2].starts_at_start?) ?
              SINGLE_BASE_REVCOM[seq[read.start_coord-1]] :
              seq[read.start_coord+middle_nodes_length]
            left_comparison_base = nodes[0].starts_at_start? ?
              first_node.ends_of_kmers_of_twin_node[0] :
              first_node.ends_of_kmers_of_node[0]
            if left_base != left_comparison_base
              log.debug "left comparison base mismatch, this is not a validating read" if log.debug?
              next
            end

            # right base, overlapping the last node
            last_node = nodes[-1].node
            right_base = !(read.direction ^ nodes[-2].starts_at_start?) ?
              seq[read.start_coord+middle_nodes_length] :
              SINGLE_BASE_REVCOM[seq[read.start_coord-1]]
            right_comparison_base = nodes[-1].starts_at_start? ?
              last_node.ends_of_kmers_of_node[0] :
              last_node.ends_of_kmers_of_twin_node[0]
            if right_base != right_comparison_base
              log.debug "right comparison base mismatch, this is not a validating read" if log.debug?
              next
            end

            log.debug "Read validates path"
            return true #gauntlet passed, this is a confirmatory read, and so the path is validated.
          end
        end
        return false #no candidate reads pass
      end

      def find_paths_from_problems(problems, recoherence_kmer)
        stack = DS::Stack.new

        to_return = Bio::AssemblyGraphAlgorithms::TrailSet.new
        to_return.circular_paths_detected = false
        to_return.trails = []


        # if there is no solutions to the overall problem then there is no solution at all
        if problems.terminal_node_keys.nil? or problems.terminal_node_keys.empty?
          return to_return
        end

        # push all solutions to the "ending in the final node" solutions to the stack
        problems.terminal_node_keys.each do |key|
          overall_solution = problems[key]
          stack.push [
            overall_solution.known_paths[0].to_a,
            []
          ]
        end
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
              paths_to_last = problems[array_trail_to_settable(first_half, recoherence_kmer)].known_paths
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
              second_to_last = path[-2]
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

      # Like a Hash, but also contains a list of keys that end in the
      # terminal node
      class ProblemSet < Hash
        # Array of keys to this hash that end in the terminal onode
        attr_accessor :terminal_node_keys
      end
    end
  end
end

require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder
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
  #
  # Options:
  # * max_gapfill_paths: the maxmimum number of paths to return. If this maximum is exceeded, an empty solution set is returned
  def find_all_connections_between_two_nodes(graph, initial_path, terminal_oriented_node,
    leash_length, recoherence_kmer, sequence_hash, options={})

    problems = find_all_problems(graph, initial_path, terminal_oriented_node, leash_length, recoherence_kmer, sequence_hash, options)

    paths = find_paths_from_problems(problems, recoherence_kmer, options)
    return paths
  end

  # Options:
  #
  # :max_explore_nodes: only explore this many nodes, not further.
  def find_all_problems(graph, initial_path, terminal_node, leash_length, recoherence_kmer, sequence_hash, options={})
    # setup dynamic programming cache
    problems = ProblemSet.new

    # setup stack to keep track of initial nodes
    pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
    pqueue.enqueue initial_path.copy, 0

    push_next_neighbours = lambda do |current_path|
      next_nodes = current_path.neighbours_of_last_node(graph)
      log.debug "Pushing #{next_nodes.length} new neighbours of #{current_path.last}" if log.debug?
      #TODO: not neccessary to copy all paths, can just continue one of them
      next_nodes.each do |n|
        log.debug "Pushing neighbour to stack: #{n}" if log.debug?
        path = current_path.copy
        path.add_oriented_node n
        pqueue.enqueue path, path.length_in_bp
      end
    end

    current_oriented_node_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    last_number_of_problems_observed_checkpoint = 0

    while current_path = pqueue.dequeue
      path_length = current_path.length_in_bp
      log.debug "considering #{current_path}, path length: #{path_length}" if log.debug?

      # Have we solved this before? If so, add this path to that solved problem.
      set_key = path_to_settable current_path, recoherence_kmer
      log.debug "Set key is #{set_key}" if log.debug?

      # Unless the path validates, forget it.
      if recoherence_kmer.nil?
        # Continue, assume that it validates if there is no recoherence_kmer
      elsif !validate_last_node_of_path_by_recoherence(current_path, recoherence_kmer, sequence_hash)
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

        num_done = problems.length
        if num_done > 0 and num_done % 512 == 0
          log.info "So far worked with #{num_done} head node sets, up to distance #{path_length}" if log.info?
        end
        if options[:max_explore_nodes] and num_done > options[:max_explore_nodes]
          log.warn "Explored too many nodes (#{num_done}), giving up.."
          problems = ProblemSet.new
          break
        end

        # explore the forward neighbours
        push_next_neighbours.call current_path
      end
      log.debug "Priority queue size: #{pqueue.length}" if log.debug?
    end

    return problems
  end

  def path_to_settable(path, recoherence_kmer)
    log.debug "Making settable a path: #{path}" if log.debug?
    return array_trail_to_settable(path.trail, recoherence_kmer)
  end

  def array_trail_to_settable(trail, recoherence_kmer)
    return trail.last.to_settable if recoherence_kmer.nil?

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
  def validate_last_node_of_path_by_recoherence(path, recoherence_kmer, sequence_hash, min_concurring_reads=1)
    #not possible to fail on a 1 or 2 node path, by debruijn graph definition.
    #TODO: that ain't true! If one of the two nodes is sufficiently long, reads may not agree.
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
      # break if the recoherence_kmer doesn't cover
      break if length_of_nodes.call(collected_nodes) + 1 >= recoherence_kmer
    end
    log.debug "validate: Collected nodes: #{collected_nodes}" if log.debug?
    if collected_nodes.length < 2
      log.debug "Only #{collected_nodes.length+1} nodes being tested for validation, so returning validated" if log.debug?
      return true
    end

    # There should be at least 1 read that spans the collected nodes and the last node
    # The trail validates if the above statement is true.
    #TODO: there's a possible 'bug' here in that there's guarantee that the read is overlays the
    # nodes in a consecutive and gapless manner. But I suspect that is unlikely to be a problem in practice.
    final_node = path.trail[-1].node
    possible_reads = final_node.short_reads.collect{|nr| nr.read_id}
    log.debug "validate starting from #{final_node.node_id}: Initial short reads: #{possible_reads.join(',') }" if log.debug?
    collected_nodes.each do |node|
      log.debug "Validating node #{node}" if log.debug?
      current_set = Set.new node.node.short_reads.collect{|nr| nr.read_id}
      possible_reads.select! do |r|
        current_set.include? r
      end
      if possible_reads.length < min_concurring_reads
        log.debug "First line validation failed, now detecting sub-kmer sequence overlap" if log.debug?
        trail_to_validate = path.trail[i+1..-1]
        return sub_kmer_sequence_overlap?(trail_to_validate, sequence_hash, min_concurring_reads)
      end
    end
    log.debug "Found #{possible_reads.length} reads that concurred with validation e.g. #{possible_reads[0]}" if log.debug?
    return true
  end

  # Is there overlap across the given nodes, even if the overlap
  # does not include an entire kmer?
  # nodes: an OrientedNodeTrail. To validate, there must be at least 1 read that spans all of these nodes
  # sequence_hash: Bio::Velvet::Sequence object with the sequences from the reads in the nodes
  def sub_kmer_sequence_overlap?(nodes, sequence_hash, min_concurring_reads=1)
    raise if nodes.length < 3 #should not get here - this is taken care of above
    log.debug "validating by sub-kmer sequence overlap with min #{min_concurring_reads}: #{nodes}" if log.debug?

    # Only reads that are in the second last node are possible, by de-bruijn graph definition.
    candidate_noded_reads = nodes[-2].node.short_reads
    middle_nodes_length = nodes[1..-2].reduce(0){|sum, n| sum += n.node.length}+
      +nodes[0].node.parent_graph.hash_length-1
    log.debug "Found middle nodes length #{middle_nodes_length}" if log.debug?

    num_confirming_reads = 0

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
        num_confirming_reads += 1
        if num_confirming_reads >= min_concurring_reads
          return true #gauntlet passed, this is enough confirmatory reads, and so the path is validated.
        end
      end
    end
    return false #no candidate reads pass
  end

  def find_paths_from_problems(problems, recoherence_kmer, options={})
    max_num_paths = options[:max_gapfill_paths]
    max_num_paths ||= 2196
    max_cycles = options[:max_cycles] || 1

    # Separate stacks for valid paths and paths which exceed the maximum allowed
    # cycle count.
    # Each backtrack spawns a set of new paths, which are cycle counted. If any cycle
    # is repeated more than max_cycles, the new path is pushed to the max_cycle_stack,
    # otherwise the path is pushed to the main stack. Main stack paths are prioritised.
    # The parachute can kill the search once the main stack exceeds max_gapfill_paths,
    # since all paths on it are valid.
    # The max_cycle_stack paths must be tracked until cycle repeats in second_part exceed
    # max_cycles, as they can spawn valid paths with backtracking.
    stack = DS::Stack.new
    max_cycle_stack = DS::Stack.new
    counter = CycleCounter.new(max_cycles)
    to_return = Bio::AssemblyGraphAlgorithms::TrailSet.new

    # if there is no solutions to the overall problem then there is no solution at all
    if problems.terminal_node_keys.nil? or problems.terminal_node_keys.empty?
      to_return.trails = []
      return to_return
    end

    # push all solutions to the "ending in the final node" solutions to the stack
    problems.terminal_node_keys.each do |key|
      overall_solution = problems[key]
      first_part = overall_solution.known_paths[0].to_a
      stack.push [
        first_part,
        [],
        ]
    end

    all_paths_hash = {}
    while path_parts = stack.pop || max_cycle_stack.pop
      log.debug path_parts.collect{|half| half.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') if log.debug?
      first_part = path_parts[0]
      second_part = path_parts[1]

      if first_part.length == 0
        # If we've tracked all the way to the beginning,
        # then there's no need to track further

        # add this solution if required
        # I've had some trouble getting the Ruby Set to work here, but this is effectively the same thing.
        log.debug "Found solution: #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
        key = second_part.hash
        all_paths_hash[key] ||= second_part.flatten
      else
        last = first_part.last
        if second_part.include?(last)
          log.debug "Cycle at node #{last.node_id} detected in previous path #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
          to_return.circular_paths_detected = true unless to_return.circular_paths_detected
          if max_cycles == 0  or max_cycles < counter.path_cycle_count([last, second_part].flatten)
            log.debug "Not finishing cyclic path with too many repeated cycles." if log.debug?
            next
          end
        #else
          #cycle_count = 0
        end
        paths_to_last = problems[array_trail_to_settable(first_part, recoherence_kmer)].known_paths
        paths_to_last.each do |path|
          to_push = [path[0...(path.length-1)],[last,second_part].flatten]
          if max_cycles < counter.path_cycle_count(to_push.flatten)
            log.debug "Pushing #{to_push.collect{|part| part.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') } to secondary stack" if log.debug?
            max_cycle_stack.push to_push
          else
            log.debug "Pushing #{to_push.collect{|part| part.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') } to main stack" if log.debug?
            stack.push to_push
          end
        end
      end

      # max_num_paths parachute
      if !max_num_paths.nil? and (stack.size + all_paths_hash.length) > max_num_paths
        log.info "Exceeded the maximum number of allowable paths in this gapfill" if log.info?
        to_return.max_path_limit_exceeded = true
        all_paths_hash = {}
        break
      end

    end

    to_return.trails = all_paths_hash.values
    return to_return
  end

  class CycleCounter
    include Bio::FinishM::Logging

    def initialize(max_cycles, options = {})
      @max_cycles = max_cycles
      @path_cache = Hash.new # Cache max_cycles for previously seen paths
      @forward = options[:forward] || false # By default builds hash moving backwards from end of path. This flag will reverse path direction and build hash moving forwards.
    end


    # Iterate through unique nodes of path and find maximal cycle counts
    def path_cycle_count(path)
      log.debug "Finding cycles in path #{path.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
      first_part = []
      second_part = path
      keys = []
      count = nil
      reached_max_cycles = false

      secondpath = second_part.reverse if @forward

      # Iterate along path and look for the remaining path in cache. Remember the iterated
      # path and the remaining path. Stop if a cache count is found, else use zero.
      while count.nil? and !second_part.empty?
        key = second_part.collect{|onode| onode.to_settable}.flatten

        # Check if path value is cached
        if @path_cache.has_key? key
          count = @path_cache[key]
          #log.debug "Found cached count #{count} for path #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
          break
        else
          first_part = [first_part, second_part.first].flatten
          second_part = second_part[1..-1]
        end
      end

      if second_part.empty?
        #log.debug "Reached end of path without finding cached count." if log.debug?
        count = 0
      end

      # The max cycle count for a path is the largest of:
      #   I. Cycle count for initial node of path in remaining path (without initial
      #      node), or
      #   II. Max cycle count of remaining path.

      # We then iterate back through the iterated path. If count does not exceed
      # max_cycles, we count cycles for each node in the remaining path, then
      # backtrack by moving the node to the remaining path set. We record the count
      # for each remaining path
      while !first_part.empty?

        node = first_part.last
        if !reached_max_cycles
          #log.debug "Next node is #{node.node.node_id}." if log.debug?
          node_count = path_cycle_count_for_node(node, second_part, @max_cycles)
          count = [count, node_count].max
          reached_max_cycles = count > @max_cycles
        end

        second_part = [node, second_part].flatten
        first_part = first_part[0...-1]

        key = second_part.collect{|onode| onode.to_settable}.flatten
        @path_cache[key] = count
        #log.debug "Caching cycle count #{count} for path #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
      end
      if reached_max_cycles and log.debug?
        log.debug "Most repeated cycle in path occured #{count} or more times."
      elsif log.debug?
        log.debug "Most repeated cycle in path occured #{count} times."
      end
      return count
    end

    # For an initial node, find and count unique 'simple'cycles in a path that begin at the initial
    # node, up to a max_cycles. Return count for the maximally repeated cycle if less than max_cycles,
    # or max_cycles.
    def path_cycle_count_for_node(node, path, max_cycles=1)
      #log.debug "Finding all simple cycles for node #{node.node_id} in path #{path.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
      remaining = path
      cycles = Hash.new

      remaining = remaining.reverse if @forward

      while remaining.include?(node)
        position = remaining.index(node)
        cycle = remaining[0..position]
        remaining = remaining[(position+1)..-1]
        #log.debug "Found cycle: #{cycle.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?

        set_key = cycle.collect{|onode| onode.to_settable}.flatten
        cycles[set_key] ||= 0
        cycles[set_key] += 1
        #log.debug "Found repeat #{cycles[set_key]}." if log.debug?

        if cycles[set_key] > max_cycles
          #log.debug "Max cycles #{max_cycles} exceeded." if log.debug?
          return cycles[set_key]
        end
      end
      if cycles.empty?
        max_counts = 0
      else
        max_counts = cycles.values.max
      end
      #log.debug "Most cycles found #{max_counts}." if log.debug?
      return max_counts
    end
  end

  class DynamicProgrammingProblem
    attr_accessor :min_distance, :known_paths

    def initialize
      @known_paths = []
    end
  end

  # Like a Hash, but also contains a list of keys that end in the
  # terminal node
  class ProblemSet < Hash
    # Array of keys to this hash that end in the terminal onode
    attr_accessor :terminal_node_keys
  end
end



require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AllOrfsFinder
      include Bio::FinishM::Logging

      CODON_LENGTH = 3
      START_CODONS = Set.new ['ATG']
      STOP_CODONS = Set.new ['TAG', 'TAA', 'TGA']

      # Search for open reading frames in a graph, in all the paths between
      # an initial and a terminal node
      def find_foward_orfs_in_connection(graph, initial_path, terminal_node,
        minimum_orf_length=100, leash_length=nil, max_gapfill_paths=nil, max_cycles=nil, max_explore_nodes=nil)

        problems = find_all_problems(graph, initial_path, {
          :terminal_nodes => terminal_nodes,
          :leash_length => leash_length
          })

        orf_trails = find_orfs_from_problems(problems, {
          :minimum_orf_length => minimum_orf_length,
          :max_gapfill_paths => max_gapfill_paths,
          :max_cycles => max_cycles,
          :max_explore_nodes => max_explore_nodes,
          })
      end


      def find_all_problems(graph, initial_path, options={})
        problems = SingleCoherentPathsBetweenNodesFinder::ProblemSet.new
        prob_finder = SingleCoherentPathsBetweenNodesFinder::ProblemTrailFinder.new(graph, initial_path)

        while current_path = prob_finder.dequeue
          path_length = current_path.length_in_bp
          log.debug "considering #{current_path}, path length: #{path_length}" if log.debug?


          terminal = prob_finder.last_is_terminal?(current_path) || (options[:terminal_nodes] && options[:terminal_nodes].include?(current_path.last))

          unless terminal
            if !options[:leash_length].nil? and path_length > options[:leash_length]
              # we are past the leash length, give up
              log.debug "Past leash length, giving up" if log.debug?
              next
            elsif options[:max_explore_nodes] and num_done > options[:max_explore_nodes]
              log.debug "Explored enough nodes (#{num_done}), giving up" if log.debug?
              next
            end
          end

          set_key = path_to_settable(current_path)


          if problems.has_key? set_key
            log.debug "Already seen this problem" if log.debug?
            prob = problems[set_key]
            prob.known_paths.push current_path

            # I don't think this should happen?
            raise "programming error" if path_length < prob.min_distance
            next
          end

          log.debug "New dynamic problem being solved" if log.debug?
          # new problem being solved here
          problem = SingleCoherentPathsBetweenNodesFinder::DynamicProgrammingProblem.new
          problem.min_distance = path_length
          problem.known_paths.push current_path.copy
          problems[set_key] = problem

          if terminal
            log.debug "last is terminal" if log.debug?

            problems.terminal_node_keys ||= Set.new
            problems.terminal_node_keys << set_key
            next
          end


          num_done = problems.length
          if num_done > 0 and num_done % 512 == 0
            log.info "So far worked with #{num_done} head node sets, up to distance #{path_length}" if log.info?
          end

          # explore the forward neighbours
          prob_finder.push_next_neighbours current_path
          log.debug "Priority queue size: #{prob_finder.length}" if log.debug?
        end

        return problems
      end

      def path_to_settable(path)
        return SingleCoherentPathsBetweenNodesFinder.new.path_to_settable(path, path.last.node.length_alone + CODON_LENGTH - 1)
      end


      def find_orfs_from_problems(problems, options={})
        max_num_paths = options[:max_gapfill_paths]
        max_num_paths ||= 2196
        max_cycles = options[:max_cycles] || 1
        min_orf_length = options[:min_orf_length] || 0

        counter = SingleCoherentPathsBetweenNodesFinder::CycleCounter.new(max_cycles)
        decide_stack = lambda do |to_push|
          part_nodes = [to_push[0].trail, to_push[1].otrail ? to_push[1].otrail.trail : []]
          if max_cycles < counter.path_cycle_count(part_nodes.flatten)
            log.debug "Pushing #{part_nodes.collect{|part| part.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') } to secondary stack" if log.debug?
            return true
          else
            log.debug "Pushing #{part_nodes.collect{|part| part.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') } to main stack" if log.debug?
            return false
          end
        end

        stack = SingleCoherentPathsBetweenNodesFinder::DualStack.new &decide_stack
        to_return = Bio::AssemblyGraphAlgorithms::TrailSet.new

        # if there is no solutions to the overall problem then there is no solution at all
        if problems.terminal_node_keys.nil? or problems.terminal_node_keys.empty?
          to_return.trails = []
          return to_return
        end

        # push all "ending in the final node" solutions to the stack
        problems.terminal_node_keys.each do |key|
          overall_solution = problems[key]
          first_part = overall_solution.known_paths[0].copy
          second_part = ORFsTracingTrail.new
          second_part.otrail = Bio::Velvet::Graph::OrientedNodeTrail.new
          stack.push [first_part, second_part]
        end

        all_paths_hash = {}
        while path_parts = stack.pop
          first_part = path_parts[0]
          second_part = path_parts[1]
          log.debug "#{first_part.to_shorthand} and #{second_part.otrail.to_shorthand}" if log.debug?

          # If an unclosed stop has been previously found, check for positions of stops
          current_fwd_orfs = second_part.fwd_orfs_result
          current_twin_orfs = second_part.twin_orfs_result
          if (current_fwd_orfs and not current_fwd_orfs.final_stop_positions.empty?) or (current_twin_orfs and not current_twin_orfs.final_stop_positions.empty?)
            log.debug "Searching first node of second part for codons" if log.debug
            fwd_starts, twin_stops = search_for_codons(second_part.otrail, false) # search from first

            # forward direction
            if current_fwd_orfs and not current_fwd_orfs.final_stop_positions.empty? and not fwd_starts.empty?
              log.debug "Attempt to pair start codons in forward direction at #{fwd_starts.join(',')}" if log.debug?
              fwd_result = orfs_from_start_stop_indices(fwd_starts, current_fwd_orfs.final_stop_positions, min_orf_length)
              log.debug "Found pairs #{fwd_result.start_stop_pairs.collect{|pair| pair.join(',')}.join('],[')}" if log.debug?

              # collect previous start-stop pairs
              fwd_result.start_stop_pairs.concat current_fwd_orfs.start_stop_pairs
              second_part.fwd_orfs_result = fwd_result
              log.debug "Remaining forward stops: #{second_part.fwd_orfs_result.final_stop_positions}" if log.debug?
            end

            # reverse direction
            if current_twin_orfs and not current_twin_orfs.final_stop_positions.empty? and not twin_stops.empty?
              log.debug "Attempt to pair stop codons in reverse direction at #{twin_stops.join(',')}" if log.debug?
              # twin stop positons are relative to start of first path twin node
              # add length of rest of path to get position relative to start of last path twin node
              length_of_rest_of_path = second_part.otrail.length_in_bos_within_path - second_part.otrail[0].node.length_alone
              twin_stops = twin_stops.collect{|pos| pos+length_of_rest_of_path}

              twin_result = orfs_from_start_stop_indices(current_twin_orfs.final_stop_positions, twin_stops, min_orf_length)
              log.debug "Found pairs #{twin_result.start_stop_pairs.collect{|pair| pair.join(',')}.join('],[')}" if log.debug?

              # collect previous start-stop pairs
              twin_result.start_stop_pairs.concat current_twin_orfs.start_stop_pairs
              second_part.twin_orfs_result = twin_result
              log.debug "Remaining twin starts: #{second_part.twin_orfs_result.final_stop_positions}" if log.debug?
            end
          end

          if first_part.length == 0
            # If we've tracked all the way to the beginning, then there's no need to track further
            key = second_part.otrail.trail.hash
            all_paths_hash[key] ||= second_part
            next
          end

          last = first_part.last
          if second_part.otrail.trail.include? last
            log.debug "Cycle at node #{last.node_id} detected in previous path #{second_part.collect{|onode| onode.node_id}.join(',')}." if log.debug?
            to_return.circular_paths_detected = true
            if max_cycles == 0 or max_cycles < counter.path_cycle_count([last, second_part.otrail.trail].flatten)
              log.debug "Not finishing cyclic path with too many repeated cycles." if log.debug?
              next
            end
          end

          paths_to_last = problems[path_to_settable(first_part)].known_paths
          paths_to_last.each do |path|
            new_second_part = ORFsTracingTrail.new
            new_second_part.otrail = second_part.otrail.copy
            new_second_part.otrail.trail.unshift last

            if new_second_part.fwd_orfs_result
              new_second_part.fwd_orfs_result.offset last.node.length_alone
            end

            new_first_part = path.copy
            new_first_part.remove_last_node

            stack.push [new_first_part, new_second_part]
          end

          # max_num_paths parachute
          # The parachute can kill the search once the main stack exceeds max_gapfill_paths,
          # since all paths on it are valid.
          if !max_num_paths.nil? and (stack.sizes[0] + all_paths_hash.length) > max_num_paths
            log.info "Exceeded the maximum number of allowable paths in this gapfill" if log.info?
            to_return.max_path_limit_exceeded = true
            all_paths_hash = {}
            break
          end
        end

        to_return.trails = all_paths_hash.values
        return to_return
      end

      # Returns:
      # in 'forward' mode (from_end=false)
      #   array of positions relative to start of first node of ends of start codons
      #   array of positions relative to start of first node twin of end of stop codons
      # in 'backward' mode (from_end=true)
      #   array of positions relative to start of last node of ends of stop codons
      #   array of positions relative to start of last node twin of ends of start codons
      def search_for_codons(otrail, from_end=false)
        onode = from_end ? otrail.last : otrail[0]

        # search within first / last node
        fwd_nodes_sequence, twin_nodes_sequence = get_sequences onode
        if from_end
          log.debug "Looking for stops in #{fwd_nodes_sequence}" if log.debug?
          fwd_stops_within_last = word_search(fwd_nodes_sequence, STOP_CODONS, CODON_LENGTH)
          log.debug "Found stops #{fwd_stops_within_last.keys.join(',')} at positions #{fwd_stops_within_last.values.flatten.join(',')}" if log.debug?

          log.debug "Looking for starts in #{twin_nodes_sequence}" if log.debug?
          twin_starts_within_last = word_search(twin_nodes_sequence, START_CODONS, CODON_LENGTH)
          log.debug "Found starts #{twin_starts_within_last.keys.join(',')} at positinos #{twin_starts_within_last.values.flatten.join(',')}" if log.debug?
        else
          log.debug "Looking for starts in #{fwd_nodes_sequence}" if log.debug?
          fwd_starts_within_first = word_search(fwd_nodes_sequence, START_CODONS, CODON_LENGTH)
          log.debug "Found starts #{fwd_starts_within_last.keys.join(',')} at positions #{fwd_starts_within_last.values.flatten.join(',')}" if log.debug?

          log.debug "Looking for stops in #{twin_nodes_sequence}" if log.debug?
          twin_stops_within_first = word_search(twin_nodes_sequence, STOP_CODONS, CODON_LENGTH)
          log.debug "Found stops #{twin_stops_within_last.keys.join(',')} at positinos #{twin_stops_within_last.values.flatten.join(',')}" if log.debug?
        end

        # extend search along trail
        fwd_overlap_sequence, twin_overlap_sequence = get_overlap_sequences(otrail, CODON_LENGTH, from_end)
        if from_end
          log.debug "Looking for stops in #{fwd_overlap_sequence}" if log.debug?
          fwd_stops_in_overlap = word_search(fwd_overlap_sequence, STOP_CODONS, CODON_LENGTH)
          log.debug "Found stops #{fwd_stops_in_overlap.keys.join(',')} at positions #{fwd_stops_in_overlap.values.flatten.join(',')}" if log.debug?

          log.debug "Looking for starts in #{twin_overlap_sequence}" if log.debug?
          twin_starts_in_overlap = word_search(twin_overlap_sequence, START_CODONS, CODON_LENGTH)
          log.debug "Found starts #{twin_starts_in_overlap.keys.join(',')} at positions #{twin_starts_in_overlap.values.flatten.join(',')}" if log.debug?
        else
          log.debug "Looking for starts in #{fwd_overlap_sequence}" if log.debug?
          fwd_starts_in_overlap = word_search(fwd_overlap_sequence, START_CODONS, CODON_LENGTH)
          log.debug "Found starts #{fwd_starts_in_overlap.keys.join(',')} at positions #{fwd_starts_in_overlap.values.flatten.join(',')}" if log.debug?

          log.debug "Looking for stops in #{twin_overlap_sequence}" if log.debug?
          twin_stops_in_overlap = word_search(twin_overlap_sequence, STOP_CODONS, CODON_LENGTH)
          log.debug "Found stops #{twin_stops_in_overlap.keys.join(',')} at positions #{twin_stops_in_overlap.values.flatten.join(',')}" if log.debug?
        end

        # combine
        if from_end
          fwd_stops = [fwd_stops_in_overlap.values,fwd_stops_within_last.values].flatten
          twin_starts = [twin_starts_in_overlap.values,twin_starts_within_last.values].flatten
          return fwd_stops, twin_starts
        else
          # offset positions in overlap to be relative to start of node / twin node
          fwd_starts = [fwd_starts_within_first.values,fwd_starts_in_overlap.values.collect{|pos| pos+onode.length_alone}].flatten
          twin_stops = [twin_stops_within_first.values,twin_stops_in_overlap.values.collect{|pos| pos+onode.length_alone}].flatten
          # move positions to front of words
          #fwd_starts = fwd_starts.collect{|pos| pos-CODON_LENGTH}
          #twin_stops = twin_stops.collect{|pos| pos-CODON_LENGTH}
          return fwd_starts, twin_stops
        end
      end

      def get_overlap_sequences(otrail, size, from_end=false)
        return '' if otrail.trail.empty?
        trail = otrail.trail
        if from_end # reverse as new sequence is taken from front of trail
          trail = trail.reverse
        end
        twin_nodes_sequence = ''
        fwd_nodes_sequence = ''

        index = 0
        onode = trail[index]

        start_length = onode.node.length_alone
        extension_length = -start_length

        while extension_length < (size - 1) and index < trail.length
          #log.debug "Extended #{extension_length} / #{size} bps and #{index+1} / #{otrail.length} nodes" if log.debug?
          extend_fwd_nodes_sequence, extend_twin_nodes_sequence = get_sequences(onode)
          if from_end
            twin_nodes_sequence += extend_twin_nodes_sequence
            fwd_nodes_sequence = extend_fwd_nodes_sequence + fwd_nodes_sequence
          else
            twin_nodes_sequence = extend_twin_nodes_sequence + twin_nodes_sequence
            fwd_nodes_sequence += extend_fwd_nodes_sequence
          end

          extension_length += onode.node.length_alone
          index += 1
          onode = trail[index]
        end

        #log.debug "Found forward and twin sequences #{fwd_nodes_sequence} and #{twin_nodes_sequence} before trimming" if log.debug?

        trim_start = start_length - (size - 1)
        trim_start = 0 if trim_start < 0
        trim_end = extension_length - (size - 1)
        trim_end = 0 if trim_end < 0
        #log.debug "Trimming first #{trim_start} and last #{trim_end} positions for output" if log.debug?
        if from_end
          return fwd_nodes_sequence[trim_end..-(trim_start+1)], twin_nodes_sequence[trim_start..-(trim_end+1)]
        else
          return fwd_nodes_sequence[trim_start..-(trim_end+1)], twin_nodes_sequence[trim_end..-(trim_start+1)]
        end
      end

      def get_sequences(onode)
        if onode.starts_at_start?
          twin_nodes_sequence = onode.node.ends_of_kmers_of_twin_node
          fwd_nodes_sequence = onode.node.ends_of_kmers_of_node
        else
          twin_nodes_sequence = onode.node.ends_of_kmers_of_node
          fwd_nodes_sequence = onode.node.ends_of_kmers_of_twin_node
        end
        return fwd_nodes_sequence, twin_nodes_sequence
      end

      def word_search(sequence, words, size)
        position = size
        inds = {}

        while position <= sequence.length
          word = sequence[position-size...position]
          if words.include? word
            inds[word] ||=[]
            inds[word].push position
          end
          position += 1
        end

        return inds
      end


      # Given an Array of start positions and stop positions, return
      # start,stop base position pairs (not inclusive of the stop codon's bases)
      # that are ORFs with a given minimum orf length (length measured in nucleotides).
      # The returned object is an instance of ORFsResult.
      def orfs_from_start_stop_indices(start_positions, stop_positions, minimum_orf_length)
        # Split up the start and stop positions into 3 frames
        frame_starts = [[],[],[]]
        frame_stops = [[],[],[]]
        start_positions.each do |pos|
          frame_starts[pos % 3].push pos
        end
        stop_positions.each do |pos|
          frame_stops[pos % 3].push pos
        end

        # For each frame
        to_return = ORFsResult.new
        (0..2).each do |frame|
          frame_pairs = []

          # Sort arrays in descending order because Array#pop removes from the end of the array
          starts = frame_starts[frame].sort{|a,b| b<=>a}
          stops = frame_stops[frame].sort{|a,b| b<=>a}

          current_start = starts.pop
          current_stop = stops.pop
          stop_before_first_start = true

          while true
            # No more stops, done
            break if current_stop.nil?
            if stop_before_first_start and (current_start.nil? or current_stop <= current_start)
              # Found stop codon before the first start codon
              to_return.final_stop_positions.push current_stop
              stop_before_first_start = false
            end
            stop_before_first_start = false

            # If there are no more start codons, we have to do no more
            if current_start.nil?
              break
            elsif current_stop <= current_start
              current_stop = stops.pop
              next
            else
              # This stop codon stops the current reading frame.
              if current_stop-current_start >= minimum_orf_length
                # Found a legit ORF
                to_return.start_stop_pairs.push [current_start, current_stop]
              end

              # Whether or not last ORF was long enough, search for the next start codon
              while current_start and current_start < current_stop
                current_start = starts.pop
              end
            end
          end
        end

        return to_return
      end

      class ORFsTracingTrail
        attr_accessor :otrail, :fwd_orfs_result, :twin_orfs_result
        include Enumerable

        def each(&block)
          unless @otrail.nil?
            @otrail.each(&block)
          end
        end
      end

      class ORFsResult
        attr_accessor :start_stop_pairs, :final_stop_positions

        def initialize
          @start_stop_pairs = []
          @final_stop_positions = []
        end

        def offset(len)
          @start_stop_pairs.map!{|pair| pair.collect{|pos| pos+len}}
          @final_stop_positions.map!{|pos| pos+len}
        end
      end
    end
  end
end

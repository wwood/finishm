require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AllOrfsFinder
      include Bio::FinishM::Logging

      CODON_LENGTH = 3
      START_CODONS = ['ATG']
      STOP_CODONS = ['TAG', 'TAA', 'TGA']

      # Search for open reading frames in a graph, in all the paths between
      # an initial and a terminal node
      def find_orfs_in_graph(graph, minimum_orf_length=nil, initial_paths=nil,
          terminal_nodes=nil, max_gapfill_paths=nil, max_cycles=nil)

        problems = find_all_problems(graph,
          initial_paths,
          :terminal_nodes => terminal_nodes
          )

        find_orfs_from_problems(problems, {
          :min_orf_length => minimum_orf_length,
          :max_gapfill_paths => max_gapfill_paths,
          :max_cycles => max_cycles,
          })
      end


      def find_all_problems(graph, initial_paths, options={})
        problems = SingleCoherentPathsBetweenNodesFinder::ProblemSet.new
        prob_finder = AllProblemTrailsFinder.new(graph, initial_paths)

        while current_path = prob_finder.pop
          log.debug "considering #{current_path}" if log.debug?
          set_key = path_to_settable(current_path)

          if problems.has_key? set_key
            log.debug "Already seen this problem" if log.debug?
            prob = problems[set_key]
            prob.known_paths.push current_path
            next
          end

          log.debug "New dynamic problem being solved" if log.debug?
          # new problem being solved here
          problem = SingleCoherentPathsBetweenNodesFinder::DynamicProgrammingProblem.new
          problem.known_paths.push current_path.copy
          problems[set_key] = problem


          terminal = current_path.neighbours_of_last_node(graph).empty? || (options[:terminal_nodes] && options[:terminal_nodes].include?(current_path.last.node))
          if terminal
            log.debug "last is terminal" if log.debug?

            problems.terminal_node_keys ||= Set.new
            problems.terminal_node_keys << set_key
            next
          end

          # explore the forward neighbours
          prob_finder.push_next_neighbours current_path
          log.debug "Priority queue size: #{prob_finder.size}" if log.debug?
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
        min_orf_length = options[:minimum_orf_length] || 0

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

          # Look for codons
          log.debug "Searching first node of second part for codons" if log.debug?
          fwd_result, twin_result = search_for_codons(second_part.otrail) # search from start of second part

          # Forward direction
          if not fwd_result.stop_positions.empty? or not fwd_result.start_positions.empty?
            current_fwd_stops = []
            if second_part.fwd_orfs_result
              current_fwd_stops.concat second_part.fwd_orfs_result.initial_stop_positions
            end
            current_fwd_stops.concat fwd_result.stop_positions
            log.debug "Attempt to pair start codons at #{fwd_result.start_positions.join(',')}" if log.debug?
            fwd_orfs_result = orfs_from_start_stop_indices(fwd_result.start_positions, current_fwd_stops, min_orf_length)
            log.debug "Found pairs #{fwd_orfs_result.start_stop_pairs.collect{|pair| pair.join(',')}.join('],[')}" if log.debug?

            # collect previous start-stop pairs
            if second_part.fwd_orfs_result
              fwd_orfs_result.start_stop_pairs.concat second_part.fwd_orfs_result.start_stop_pairs
            end
            second_part.fwd_orfs_result = fwd_orfs_result
            log.debug "Remaining forward stops: #{second_part.fwd_orfs_result.initial_stop_positions}" if log.debug?
          end

          if not twin_result.stop_positions.empty? or not twin_result.start_positions.empty?
            current_twin_starts = []
            if second_part.twin_orfs_result
              current_twin_starts.concat second_part.twin_orfs_result.final_start_positions
            end
            current_twin_starts.concat twin_result.start_positions
            log.debug "Attempt to pair stop codons in reverse direction at #{twin_stops.join(',')}" if log.debug?
            # twin stop positons are relative to start of first path twin node
            # add length of rest of path to get position relative to start of last path twin node
            length_of_rest_of_path = second_part.otrail.length_in_bos_within_path - second_part.otrail[0].node.length_alone
            current_twin_stops = twin_result.stop_positions.collect{|pos| pos+length_of_rest_of_path}

            twin_orfs_result = orfs_from_start_stop_indices(current_twin_starts, current_twin_stops, min_orf_length)
            log.debug "Found pairs #{twin_orfs_result.start_stop_pairs.collect{|pair| pair.join(',')}.join('],[')}" if log.debug?

            # collect previous start-stop pairs
            if second_part.twin_orfs_result
              twin_orfs_result.start_stop_pairs.concat second_part.twin_orfs_result.start_stop_pairs
            end
            second_part.twin_orfs_result = twin_orfs_result
            log.debug "Remaining twin starts: #{second_part.twin_orfs_result.final_start_positions}" if log.debug?
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

            if second_part.fwd_orfs_result
              # offset positions in forward direction
              offset = last.node.length_alone
              new_fwd_orfs_result = ORFsResult.new
              new_fwd_orfs_result.start_stop_pairs = second_part.fwd_orfs_result.start_stop_pairs.collect{|pairs| pairs.collect{|pos| pos + offset}}
              new_fwd_orfs_result.initial_stop_positions = second_part.fwd_orfs_result.initial_stop_positions.collect{|pos| pos + offset}
              new_fwd_orfs_result.final_start_positions = second_part.fwd_orfs_result.final_start_positions.collect.to_a
              new_second_part.fwd_orfs_result = new_fwd_orfs_result
            end

            if second_part.twin_orfs_result
              offset ||= last.node.length_alone
              new_twin_orfs_result = ORFsResult.new
              new_twin_orfs_result.start_stop_pairs = second_part.twin_orfs_result.start_stop_pairs.collect{|pairs| pairs.collect{|pos| pos}}
              new_twin_orfs_result.initial_stop_positions = second_part.twin_orfs_result.initial_stop_positions.collect.to_a
              new_twin_orfs_result.final_start_positions = second_part.twin_orfs_result.final_start_positions.collect{|pos| pos + offset}
              new_second_part.twin_orfs_result = new_twin_orfs_result
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
      # SearchResult relative to start of first node
      # SearchResult relative to start of first twin node
      def search_for_codons(otrail)
        return SearchResult.new, SearchResult.new if otrail.trail.empty?
        onode = otrail[0]

        # search within first / last node
        fwd_nodes_sequence, twin_nodes_sequence = get_sequences onode
        words = Set.new(START_CODONS).merge(STOP_CODONS)
        #log.debug "Looking for codons #{words.to_a}" if log.debug?

        # forward
        #log.debug "Looking in #{fwd_nodes_sequence}" if log.debug?
        fwd_within_first = word_search(fwd_nodes_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{fwd_within_first.keys.join(',')} at positions #{fwd_within_first.values.join(',')} in #{fwd_nodes_sequence}" if log.debug?

        # reverse
        #log.debug "Looking in #{twin_nodes_sequence}" if log.debug?
        twin_within_first = word_search(twin_nodes_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{twin_within_first.keys.join(',')} in twin node at positions #{twin_within_first.values.join(',')} in #{fwd_nodes_sequence}" if log.debug?


        # extend search along trail
        fwd_overlap_sequence, twin_overlap_sequence = get_overlap_sequences(otrail, CODON_LENGTH)

        # forward
        #log.debug "Looking in #{fwd_overlap_sequence}" if log.debug?
        fwd_in_overlap = word_search(fwd_overlap_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{fwd_in_overlap.keys.join(',')} in twin node at positions #{fwd_in_overlap.values.join(',')} in #{fwd_overlap_sequence}" if log.debug?


        # reverse
        #log.debug "Looking for stops in #{twin_overlap_sequence}" if log.debug?
        twin_in_overlap = word_search(twin_overlap_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{twin_in_overlap.keys.join(',')} in twin node at positions #{twin_in_overlap.values.join(',')} in #{twin_overlap_sequence}" if log.debug?

        # offset positions in overlap to be relative to start of node / twin node
        offset = onode.node.length_alone
        fwd_in_overlap.each{|word, inds| fwd_in_overlap[word] = inds.collect{|pos| pos + offset}}
        twin_in_overlap.each{|word, inds| twin_in_overlap[word] = inds.collect{|pos| pos + 1 - CODON_LENGTH}}
        #log.debug "Codons in overlap positions relative to start of first node #{fwd_in_overlap.values.join(',')}" if log.debug?
        #log.debug "Codons in overlap positions relative to start of first twin node #{twin_in_overlap.values.join(',')}" if log.debug?

        # assemble result
        fwd_result = SearchResult.new
        twin_result = SearchResult.new
        START_CODONS.each do |word|
          # fwd starts
          if fwd_within_first.has_key? word
            fwd_result.start_positions.push fwd_within_first[word]
          end
          if fwd_in_overlap.has_key? word
            fwd_result.start_positions.push fwd_in_overlap[word]
          end
          # twin starts
          if twin_within_first.has_key? word
            twin_result.start_positions.push twin_within_first[word]
          end
          if twin_in_overlap.has_key? word
            twin_result.start_positions.push twin_in_overlap[word]
          end
          end
        fwd_result.start_positions.flatten!
        twin_result.start_positions.flatten!
        #log.debug "Positions of start codons #{fwd_result.start_positions.join(',')}" if log.debug?
        #log.debug "Positions of start codons in twin node #{twin_result.start_positions.join(',')}" if log.debug?

        STOP_CODONS.each do |word|
          #fwd stops
          if fwd_within_first.has_key? word
            fwd_result.stop_positions.push fwd_within_first[word]
          end
          if fwd_in_overlap.has_key? word
            fwd_result.stop_positions.push fwd_in_overlap[word]
          end
          # twin stops
          if twin_within_first.has_key? word
            twin_result.stop_positions.push twin_within_first[word]
          end
          if twin_in_overlap.has_key? word
            twin_result.stop_positions.push twin_in_overlap[word]
          end
        end
        fwd_result.stop_positions.flatten!
        twin_result.stop_positions.flatten!
        #log.debug "Positions of stop codons #{fwd_result.stop_positions.join(',')}" if log.debug?
        #log.debug "Positions of stop codons in twin node #{twin_result.stop_positions.join(',')}" if log.debug?


        return fwd_result, twin_result
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
          before_first_start = true

          while true
            if current_stop
              if current_start.nil? or current_stop <= current_start
                if before_first_start
                  # Found stop codon before the first start codon
                  to_return.initial_stop_positions.push current_stop
                  before_first_start = false
                end

                # If there are no more start codons, we have to do no more
                if current_start.nil?
                  break
                else
                  while current_stop and current_stop <= current_start
                    current_stop = stops.pop
                  end
                  next
                end
              else
                before_first_start = false
                # This stop codon stops the current reading frame.
                if current_stop-current_start >= minimum_orf_length
                  # Found a legit ORF
                  to_return.start_stop_pairs.push [current_start, current_stop]
                end

                # Whether or not last ORF was long enough, search for the next start codon
                while current_start and current_start < current_stop
                  current_start = starts.pop
                end
                next
              end
            elsif current_start
              # Found a start codon after last stop codon
              to_return.final_start_positions.push current_start
              break
            end
            break
          end
        end

        return to_return
      end

      # SearchResult with fields:
      # array of positions of last base of start codons
      # array of positions of last base of stop codons
      class SearchResult
        attr_accessor :start_positions, :stop_positions

        def initialize
          @start_positions = []
          @stop_positions = []
        end
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
        attr_accessor :start_stop_pairs, :final_start_positions, :initial_stop_positions

        def initialize
          @start_stop_pairs = []
          @final_start_positions = []
          @initial_stop_positions = []
        end
      end

      class AllProblemTrailsFinder
        include Bio::FinishM::Logging

        def initialize(graph, initial_paths)
          @stack = DS::Stack.new
          initial_paths.each do |path|
            @stack.push path
          end
          @graph = graph
        end

        def pop
          @stack.pop
        end

        def size
          @stack.size
        end

        def push_next_neighbours(current_path)
          next_nodes = current_path.neighbours_of_last_node(@graph)
          log.debug "Pushing #{next_nodes.length} new neighbours of #{current_path.last}" if log.debug?
          #TODO: not neccessary to copy all paths, can just continue one of them
          next_nodes.each do |n|
            log.debug "Pushing neighbour to stack: #{n}" if log.debug?
            path = current_path.copy
            path.add_oriented_node n
            @stack.push path
          end
        end
      end
    end
  end
end

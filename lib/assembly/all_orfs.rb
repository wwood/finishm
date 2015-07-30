require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AllOrfsFinder
      include Bio::FinishM::Logging

      CODON_LENGTH = 3
      START_CODONS = ['ATG']
      STOP_CODONS = ['TAG', 'TAA', 'TGA']
      CODONS = {
        'A' => ['GCT', 'GCC', 'GCA', 'GCG'],
        'R' => ['CGT', 'CGC', 'CGA','CGG', 'AGA', 'AGG'],
        'N' => ['AAT', 'AAC'],
        'D' => ['GAT', 'GAC'],
        'C' => ['TGT', 'TGC'],
        'Q' => ['CAA', 'CAG'],
        'E' => ['GAA', 'GAG'],
        'G' => ['GGT', 'GGC', 'GGA', 'GGG'],
        'H' => ['CAT', 'CAC'],
        'I' => ['ATT', 'ATC', 'ATA'],
        'L' => ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'K' => ['AAA', 'AAG'],
        'M' => ['ATG'],
        'F' => ['TTT', 'TTC'],
        'P' => ['CCT', 'CCC', 'CCA', 'CCG'],
        'S' => ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T' => ['ACT', 'ACC', 'ACA', 'ACG'],
        'W' => ['TGG'],
        'Y' => ['TAT', 'TAC'],
        'V' => ['GTT', 'GTC', 'GTA', 'GTG']
        }
      TRANSLATOR = CODONS.reduce({}) do |memo, pair|
        pair[1].each{|key| memo[key] = pair[0]}
        memo
      end

      # Search for open reading frames in a graph, in all the paths begining at a set of
      # nodes through a graph (or a subset defined by range)
      def find_orfs_in_graph(graph, initial_paths, minimum_orf_length=nil,
          range=nil, max_gapfill_paths=nil, max_cycles=nil)

        problems = find_all_problems(graph,
          initial_paths,
          :range => range
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

          neighbours = current_path.neighbours_of_last_node(graph)
          if options[:range]
            neighbours.select!{|onode| options[:range].include? onode.node_id}
          end
          if neighbours.empty?
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

          # Look for codons
          log.debug "Searching for codons in first node of second part" if log.debug?
          fwd_result, twin_result = search_for_codons(second_part.otrail) # search from start of second part

          # Forward direction
          if not fwd_result.stop_markers.empty? or not fwd_result.start_markers.empty?
            [fwd_result.stop_markers, fwd_result.start_markers].each do |markers|
              markers.each do |marker|
                marker.position_in_trail = marker.position_in_node
              end
            end
            current_fwd_stops = []
            current_fwd_starts = []
            if second_part.fwd_orfs_result
              current_fwd_stops.concat second_part.fwd_orfs_result.initial_stop_markers
              current_fwd_starts.concat second_part.fwd_orfs_result.initial_start_markers
              current_fwd_starts.concat second_part.fwd_orfs_result.final_start_markers
            end
            current_fwd_stops.concat fwd_result.stop_markers
            current_fwd_starts.concat fwd_result.start_markers
            log.debug "Attempt to pair start codons at #{current_fwd_starts.collect{|m| m.position_in_trail}.join(',')} with stop codons at #{current_fwd_stops.collect{|m| m.position_in_trail}.join(',')}" if log.debug?
            fwd_orfs_result = orfs_from_start_stop_markers(current_fwd_starts, current_fwd_stops, min_orf_length)
            log.debug "Found pairs #{fwd_orfs_result.start_stop_pairs.collect{|pair| pair.collect{|m| m.position_in_trail}.join(',')}.join('],[')}" if log.debug?

            # collect previous start-stop pairs
            if second_part.fwd_orfs_result
              fwd_orfs_result.start_stop_pairs.concat second_part.fwd_orfs_result.start_stop_pairs
            end
            second_part.fwd_orfs_result = fwd_orfs_result
            log.debug "Remaining forward stops: #{second_part.fwd_orfs_result.initial_stop_markers.collect{|m| m.position_in_trail}.join(',')}" if log.debug?
          end

          # Reverse direction
          if not twin_result.stop_markers.empty? or not twin_result.start_markers.empty?
            # twin stop positons are relative to start of first path twin node
            # add length of rest of path to get position relative to start of last path twin node
            length_of_rest_of_path = second_part.otrail.length_in_bp_within_path - second_part.otrail[0].node.length_alone
            [twin_result.stop_markers, twin_result.start_markers].each do |markers|
              markers.each do |marker|
                marker.position_in_trail = marker.position_in_node + length_of_rest_of_path
              end
            end
            current_twin_stops = []
            current_twin_starts = []
            if second_part.twin_orfs_result
              current_twin_stops.concat second_part.twin_orfs_result.initial_stop_markers
              current_twin_starts.concat second_part.twin_orfs_result.initial_start_markers
              current_twin_starts.concat second_part.twin_orfs_result.final_start_markers
            end
            current_twin_stops.concat twin_result.stop_markers
            current_twin_starts.concat twin_result.start_markers
            log.debug "Attempt to pair stop codons in reverse direction at #{current_twin_stops.collect{|m| m.position_in_trail}.join(',')} with starts at #{current_twin_starts.collect{|m| m.position_in_trail}.join(',')}" if log.debug?
            twin_orfs_result = orfs_from_start_stop_markers(current_twin_starts, current_twin_stops, min_orf_length)
            log.debug "Found pairs #{twin_orfs_result.start_stop_pairs.collect{|pair| pair.collect{|m| m.position_in_trail}.join(',')}.join('],[')}" if log.debug?

            # collect previous start-stop pairs
            if second_part.twin_orfs_result
              twin_orfs_result.start_stop_pairs.concat second_part.twin_orfs_result.start_stop_pairs
            end
            second_part.twin_orfs_result = twin_orfs_result
            log.debug "Remaining twin starts: #{second_part.twin_orfs_result.final_start_markers.collect{|m| m.position_in_trail}.join(',')}" if log.debug?
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
              copy_and_offset_marker = lambda do |marker|
                m = marker.copy
                m.position_in_trail += offset
                m
              end

              new_fwd_orfs_result = ORFsResult.new
              new_fwd_orfs_result.start_stop_pairs = second_part.fwd_orfs_result.start_stop_pairs.collect do |pairs|
                pairs.collect &copy_and_offset_marker
              end
              new_fwd_orfs_result.initial_start_markers = second_part.fwd_orfs_result.initial_start_markers.collect &copy_and_offset_marker
              new_fwd_orfs_result.initial_stop_markers = second_part.fwd_orfs_result.initial_stop_markers.collect &copy_and_offset_marker
              new_fwd_orfs_result.final_start_markers = second_part.fwd_orfs_result.final_start_markers.collect &copy_and_offset_marker
              new_second_part.fwd_orfs_result = new_fwd_orfs_result
            end

            if second_part.twin_orfs_result
              new_twin_orfs_result = ORFsResult.new
              new_twin_orfs_result.start_stop_pairs = second_part.twin_orfs_result.start_stop_pairs.collect do |pairs|
                pairs.collect{|marker| marker.copy}
              end
              new_twin_orfs_result.initial_stop_markers = second_part.twin_orfs_result.initial_stop_markers.collect{|marker| marker.copy}
              new_twin_orfs_result.final_start_markers = second_part.twin_orfs_result.final_start_markers.collect{|marker| marker.copy}
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

        make_marker = lambda do |position|
          marker = Marker.new
          marker.position_in_node = position
          marker.node = onode.node
          marker
        end

        #log.debug "Looking for codons #{words.to_a}" if log.debug?
        words = Set.new(START_CODONS).merge(STOP_CODONS)

        # search within first / last node
        fwd_nodes_sequence, twin_nodes_sequence = get_sequences onode
        #log.debug "Looking in #{fwd_nodes_sequence}" if log.debug?
        fwd_within_first = word_search(fwd_nodes_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{fwd_within_first.keys.join(',')} at positions #{fwd_within_first.values.join(',')} in #{fwd_nodes_sequence}" if log.debug?
        #log.debug "Looking in #{twin_nodes_sequence}" if log.debug?
        twin_within_first = word_search(twin_nodes_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{twin_within_first.keys.join(',')} in twin node at positions #{twin_within_first.values.join(',')} in #{fwd_nodes_sequence}" if log.debug?

        # extend search along trail
        fwd_overlap_sequence, twin_overlap_sequence = get_overlap_sequences(otrail, CODON_LENGTH)
        #log.debug "Looking in #{fwd_overlap_sequence}" if log.debug?
        fwd_in_overlap = word_search(fwd_overlap_sequence, words, CODON_LENGTH)
        #log.debug "Found codons #{fwd_in_overlap.keys.join(',')} in twin node at positions #{fwd_in_overlap.values.join(',')} in #{fwd_overlap_sequence}" if log.debug?
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

        push_mark_to_list = lambda do |list, word, positions|
          if positions.has_key? word
            list.push positions[word].collect{|pos| make_marker.call pos}
          end
        end

        fwd_positions = fwd_within_first.merge fwd_in_overlap
        twin_positions = twin_within_first.merge twin_in_overlap
        START_CODONS.each do |word|
          # fwd starts
          push_mark_to_list.call(fwd_result.start_markers, word, fwd_positions)
          # twin starts
          push_mark_to_list.call(twin_result.start_markers, word, twin_positions)
        end
        fwd_result.start_markers.flatten!
        twin_result.start_markers.flatten!
        #log.debug "Positions of start codons #{fwd_result.start_markers.join(',')}" if log.debug?
        #log.debug "Positions of start codons in twin node #{twin_result.start_markers.join(',')}" if log.debug?

        STOP_CODONS.each do |word|
          #fwd stops
          push_mark_to_list.call(fwd_result.stop_markers, word, fwd_positions)
          # twin stops
          push_mark_to_list.call(twin_result.stop_markers, word, twin_positions)
        end
        fwd_result.stop_markers.flatten!
        twin_result.stop_markers.flatten!
        #log.debug "Positions of stop codons #{fwd_result.stop_markers.join(',')}" if log.debug?
        #log.debug "Positions of stop codons in twin node #{twin_result.stop_markers.join(',')}" if log.debug?


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
      def orfs_from_start_stop_markers(start_markers, stop_markers, minimum_orf_length)
        # Split up the start and stop positions into 3 frames
        frame_starts = [[],[],[]]
        frame_stops = [[],[],[]]
        start_markers.each do |marker|
          frame_starts[marker.position_in_trail % 3].push marker
        end
        stop_markers.each do |marker|
          frame_stops[marker.position_in_trail % 3].push marker
        end

        # For each frame
        to_return = ORFsResult.new
        (0..2).each do |frame|
          frame_pairs = []

          # Sort arrays in descending order because Array#pop removes from the end of the array
          starts = frame_starts[frame].sort{|a,b| b.position_in_trail<=>a.position_in_trail}
          stops = frame_stops[frame].sort{|a,b| b.position_in_trail<=>a.position_in_trail}

          current_start = starts.pop
          current_stop = stops.pop
          if current_stop
            # Record first stop codon
            to_return.initial_stop_markers.push current_stop
          end
          if current_start and (current_stop.nil? or current_start.position_in_trail < current_stop.position_in_trail)
            # Record first start codon before any stop codons
            to_return.initial_start_markers.push current_start
          end

          while current_start and current_stop
            # Move to next start after current stop
            while current_start and current_start.position_in_trail < current_stop.position_in_trail
              current_start = starts.pop
            end

            if current_start
              # Move to next stop after current start
              while current_stop and current_stop.position_in_trail < current_start.position_in_trail
                current_stop = stops.pop
              end
            end

            if current_start and current_stop
              # This stop codon stops the current reading frame.
              if current_stop.position_in_trail - current_start.position_in_trail >= minimum_orf_length
                # Found a legit ORF
                to_return.start_stop_pairs.push [current_start, current_stop]
              end
              # Whether or not last ORF was long enough, search for the next start codon
              next
            else
              if current_start
                to_return.final_start_markers.push current_start
              end
              break
            end
          end
        end

        return to_return
      end

      def orf_sequences_from_trails(trails, min_orf_length=nil)
        to_return = []
        trails.each do |trail|
          fwd_sequence, twin_sequence = trail.otrail.sequences_within_path
          trail_length = fwd_sequence.length
          # forward / twin directions
          [
          [fwd_sequence, trail.fwd_orfs_result],
          [twin_sequence, trail.twin_orfs_result]
          ].each do |sequence_and_result|
            sequence, result = sequence_and_result
            if result
              result.start_stop_pairs.each do |pair|
                start_position = pair[0].position_in_trail - 3
                end_position = pair[1].position_in_trail

                # orf name
                last_node = nil
                onodes = trail.otrail.trail.drop_while do |onode|
                  onode.node != pair[0].node
                end.take_while do |onode|
                  next false if last_node == pair[1].node
                  last_node = onode.node
                  true
                end
                name = "(#{onodes[0].to_shorthand}:#{pair[0].position_in_node}),#{onodes[1...-1].collect{|onode| onode.to_shorthand}.join(',')},(#{onodes[-1].to_shorthand}:#{pair[1].position_in_node})"

                to_return.push [name, sequence[start_position...end_position]]
              end
              result.initial_stop_markers.each do |marker|
                end_position = marker.position_in_trail
                start_position = end_position % 3 #trim sequence to multiple of 3
                next if min_orf_length and end_position - start_position < min_orf_length

                # orf_name
                last_node = nil
                onodes = trail.otrail.trail.take_while do |onode|
                  next false if last_node == marker.node
                  last_node = onode.node
                  true
                end
                name = "#{onodes[0...-1].collect{|onode| onode.to_shorthand}.join(',')},(#{onodes[-1].to_shorthand}:#{marker.position_in_node})"

                to_return.push [name, sequence[start_position...end_position]]
              end
              result.final_start_markers.each do |marker|
                start_position = marker.position_in_trail - 3
                end_position = (trail_length - start_position) % 3
                next if min_orf_length and trail_length - end_position - start_position < min_orf_length

                # orf_name
                onodes = trail.otrail.trail.drop_while{|onode| onode.node != marker.node}
                name = "(#{onodes[0].to_shorthand}:#{marker.position_in_node}),#{onodes[1..-1].collect{|onode| onode.to_shorthand}.join(',')}"
                to_return.push [name, sequence[start_position..-1-end_position]]
              end
            end
            if result.nil? or (result.start_stop_pairs.empty? and result.final_start_markers.empty? and result.initial_stop_markers.empty?)
              (0..2).each do |frame|
                start_position = frame
                end_position = (trail_length - start_position) % 3
                next if min_orf_length and trail_length - end_position - start_position < min_orf_length

                # orf_name
                name = "#{trail.otrail.to_shorthand}"
                to_return.push [name, sequence[start_position..-1-end_position]]
              end
            end
          end
        end

        return to_return
      end

      def sequence2AA(sequence)
        remaining = sequence
        aa = ""
        while remaining.length > 0
          codon = remaining[0...3]
          log.debug "Found next codon #{codon}" if log.debug?
          if not TRANSLATOR.has_key?(codon)
            raise "Cannot translate invalid codon #{codon} in sequence #{sequence}."
          end
          log.debug "Codon translated to #{TRANSLATOR[codon]}" if log.debug?
          aa += TRANSLATOR[codon]
          remaining = remaining[3..-1]
        end
        return aa
      end

      # positions of last base of codons
      class Marker
        attr_accessor :position_in_trail, :position_in_node, :node

        def copy
          copy = Marker.new
          copy.position_in_trail = @position_in_trail
          copy.position_in_node = @position_in_node
          copy.node = @node
          return copy
        end
      end

      class SearchResult
        attr_accessor :start_markers, :stop_markers

        def initialize
          @start_markers = []
          @stop_markers = []
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
        attr_accessor :start_stop_pairs, :final_start_markers, :initial_start_markers, :initial_stop_markers

        def initialize
          @start_stop_pairs = []
          @initial_start_markers = []
          @final_start_markers = []
          @initial_stop_markers = []
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

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
        minimum_orf_length=100, leash_length=nil)
        raise "this method is not yet implemented properly. Need to write something to
* take a path and last start points of ORFs, and update them, returning ORFs
* find a way of hashing these kinds of objects so that they can be saved for dynamic programming purposes
* test, of course"

        # Setup and add initial path to the stack
        stack = DS::Stack.new
        stack.push initial_path

        while current_node = stack.pop
          # If we are circular, give up
          if raise
            # If we are beyond leash, then give up
            # Else look for new orfs given this current path
          else
            # Find all stop codons in the last node, including those that include 1 or 2 bases from the previous node, where the stop codon is across 2 (or even 3) nodes
            # Similarly, find all start codons in the last node + previous 2bp
            # TODO: maybe is somewhat wasteful to be searching for start codons in places where there is no stop codons, as they would just be internal ATGs. Therefore
            # Update the 3 ORF counters given the start and the end points
            # Proceed greedily, taking the first reading frame of the

            # record the fact that we've solved this problem (dynamic programming)

            # Go to neighbours
            # if the neighbours have already been solved, we don't need to follow those
            # trails because dynamic programming will take care of everything for us.
          end
        end
      end


      def find_all_problems(graph, initial_path, options={})
        problems = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder::ProblemSet.new
        prob_finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder::ProblemTrailFinder.new(graph, initial_path)
        path_finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

        while current_path = prob_finder.dequeue
          path_length = current_path.length_in_bp
          log.debug "considering #{current_path}, path length: #{path_length}" if log.debug?

          if prob_finder.ended? current_path  or (options[:terminal_nodes] and options[:terminal_nodes].any?{|onode| current_path.last == onode})
            log.debug "last is terminal" if log.debug?
            next
          elsif !leash_length.nil? and path_length > leash_length
            # we are past the leash length, give up
            log.debug "Past leash length, giving up" if log.debug?
            next
          elsif options[:max_explore_nodes] and num_done > options[:max_explore_nodes]
            log.debug "Explored enough nodes (#{num_done}), giving up" if log.debug?
            next
          end

          node_length_with_overlap = current_path.last.node.length_alone + CODON_LENGTH - 1
          set_key = path_finder.path_to_settable current_path, node_length_with_overlap

          if problems.has_key? set_key
            log.debug "Already seen this problem" if log.debug?
            prob = problems[set_key]
            prob.known_paths.push current_path

            # I don't think this should happen?
            raise "programming error" if path_length < prob.min_distance
            next
          end

          fwd_stops, twin_starts = search_for_words_from_last(current_path, START_CODONS, STOP_CODONS)


          log.debug "New dynamic problem being solved" if log.debug?
          # new problem being solved here
          problem = DynamicProgrammingProblem.new
          problem.min_distance = path_length
          problem.known_paths.push current_path.copy
          problems[set_key] = problem

          if not fwd_stops.empty? or not twin_starts.empty?
            problems.fwd_stops = fwd_stops unless fwd_stops.empty?
            problems.twin_starts = twin_starts unless twin_starts.empty?
            problems.terminal_node_keys ||= Set.new
            problems.terminal_node_keys << set_key
          end

          num_done = problems.length
          if num_done > 0 and num_done % 512 == 0
            log.info "So far worked with #{num_done} head node sets, up to distance #{path_length}" if log.info?
          end

          # explore the forward neighbours
          prob_finder.push_next_neighbours current_path
          log.debug "Priority queue size: #{finder.length}" if log.debug?
        end
      end


      def find_orfs_from_problems(problems, options={})
        max_num_paths = options[:max_gapfill_paths]
        max_num_paths ||= 2196
        max_cycles = options[:max_cycles] || 1

        solver = SingleCoherentPathsBetweenNodesFinder::ProblemSolver.new(problems, CODON_LENGTH, max_cycles)
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
          solver.push first_part
        end

        all_paths_hash = {}
        while paths_parts = solver.pop_parts
          log.debug path_parts.collect{|half| half.collect{|onode| onode.node.node_id}.join(',')}.join(' and ') if log.debug?
          first_part = path_parts[0]
          second_part = path_parts[1]

          if first_part.length == 0
            # If we've tracked all the way to the beginning,
            # then there's no need to track further
            next
          end

          second_part_otrail = Bio::Velvet::Graph::OrientedNodeTrail.new
          second_part_otrail.trail = second_part

          fwd_starts, twin_stops = search_for_words_from_start(second_part_otrail, START_CODONS, STOP_CODONS)

          if not fwd_starts.empty?
            # I've had some trouble getting the Ruby Set to work here, but this is effectively the same thing.
            log.debug "Found solution: #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
            next
          else
            solver.push_parts first_part, second_part
          end

          # max_num_paths parachute
          # The parachute can kill the search once the main stack exceeds max_gapfill_paths,
          # since all paths on it are valid.
          if !max_num_paths.nil? and (solver.stack_size + all_paths_hash.length) > max_num_paths
            log.info "Exceeded the maximum number of allowable paths in this gapfill" if log.info?
            to_return.max_path_limit_exceeded = true
            all_paths_hash = {}
            break
          end
        end

      end

      def search_for_words_from_last(otrail, start_words, stop_words)
        length_with_overlap = otrail.last.node.length_alone + CODON_LENGTH - 1
        fwd_sequence, twin_sequence = sequence_with_length(otrail, length_with_overlap)
        fwd_stops = stop_words.nil? ? {} : word_search(fwd_sequence, stop_words)
        twin_starts = start_words.nil? ? {} : word_search(twin_sequence.reverse, start_words)

        return fwd_stops, twin_starts
      end

      def search_for_words_from_first(otrail, start_words, stop_words)
        length_with_overlap = otrail[0].node.length_alone + CODON_LENGTH - 1
        fwd_sequence, twin_sequence = sequence_with_length(otrail, length_with_overlap, from_start=true)
        fwd_starts = start_words.nil? ? [] : word_search(fwd_sequence, start_words)
        twin_stops = stop_words.nil? ? [] : word_search(twin_sequence, stop_words)

        return fwd_starts, twin_stops
      end

      def sequence_with_length(otrail, length, from_start=false)
        return '' if otrail.trail.empty?
        if from_start
          otrail = otrail.reverse
        end

        twin_nodes_sequence = ''
        fwd_nodes_sequence = ''
        index = otrail.length - 1

        while fwd_nodes_sequence.length < length and index > 0
          onode = otrail[index]
          if onode.starts_at_start?
            twin_nodes_sequence += onode.node.ends_of_kmers_of_twin_node
            fwd_nodes_sequence = onode.node.ends_of_kmers_of_node + fwd_nodes_sequence
          else
            twin_nodes_sequence += onode.node.ends_of_kmers_of_node
            fwd_nodes_sequence = onode.node.ends_of_kmers_of_twin_node + fwd_nodes_sequence
          end
          index -= 1
        end

        return fwd_nodes_sequence[-length..-1], twin_nodes_sequence[0...length]
      end

      def word_search(padded_sequence, words)
        position = CODON_LENGTH
        inds = {}

        while position < padded_sequence.length
          # extend sequences if shorter than required

          word = padded_sequence[position-CODON_LENGTH..position]
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

          # Reverse the arrays because Array#pop removes from the end of the array
          starts = frame_starts[frame].reverse
          stops = frame_stops[frame].reverse

          current_start = starts.pop
          current_stop = nil
          if current_start.nil?
            # no start codons in the entire frame, so don't do anything
            # TODO: test what happens if there is a start but not a stop codon in this frame
          else
            while current_stop = stops.pop
              if current_stop.nil?
                # no stop codons in this frame at all. There's no ORFS detected, and the current_start is the biggest start
                to_return.final_start_positions[frame1] = current_start
                break
              elsif current_stop <= current_start
                # Skip of stop codons before the current start codon
                next
              else
                # This stop codon stops the current reading frame.
                if current_stop-current_start >= minimum_orf_length
                  # Found a legit ORF
                  to_return.start_stop_pairs[frame].push [current_start, current_stop]
                end

                # Whether or not last ORF was long enough, search for the next start codon
                while current_start = starts.pop
                  if current_start < current_stop
                    next
                  else
                    break
                  end
                end

                # If there are no more start codons, the final start position is nil, and we
                # have to do no more
                break if current_start.nil?
              end
            end

            # If there is a dangling current_start, then go from there
            to_return.final_start_positions[frame] = current_start
          end
        end

        return to_return
      end

      class ORFsResult
        attr_accessor :start_stop_pairs, :final_start_positions

        def initialize
          @start_stop_pairs = [[],[],[]]
          @final_start_positions = [nil]*3
        end
      end

      class DynamicProgrammingProblem < SingleCoherentPathsBetweenNodesFinder::DynamicProgrammingProblem
        attr_accessor :fwd_stops, :twin_starts
      end
    end
  end
end

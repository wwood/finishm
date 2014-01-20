require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AllOrfsFinder
      include Bio::FinishM::Logging

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
    end
  end
end

require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class Fluffer
      include Bio::FinishM::Logging

      class FlufferHalfResult
        attr_accessor :golden_paths, :golden_path_fragments

        def initialize
          @golden_paths = []
          @golden_path_fragments = []
        end
      end

      # Return an array of array of paths.
      def fluff(finishm_graph, leash_length, options={})
        log.debug "Fluffing part 1.."
        half_results = fluff_part1(finishm_graph, leash_length, options)

        log.debug "Found fluff half results: #{half_results}" if log.debug?
        log.debug "Fluffing part 2.."
        return fluff_part2(half_results)
      end

      def fluff_part1(finishm_graph, leash_length, options={})
        # Get a set of all probes so that they can be checked against
        log.debug "Found #{finishm_graph.probe_nodes.reject{|node| node.nil?}.length} different probes that were in the final velvet graph"
        half_results = []
        graph = finishm_graph.graph

        # For each probe in the graph
        finishm_graph.probe_nodes.each_with_index do |node, probe_index|

          # If the node is not found in the graph, then forget it
          if node.nil?
            half_results.push FlufferHalfResult.new

          else # probe was found in the graph. Start finding them paths.
            # start exploring the squid way
            golden_paths = [] #These is the full paths found
            golden_fragments = [] #paths to join up to the end
            already_visited_nodes = Set.new #Nodes that have already been explored

            golden_onodes = Set.new #Nodes that stop exploration
            terminal_nodes = Set.new
            # Add all the nodes that are being probed, because we don't want double exploration
            # i.e. want probe1 => probe2, probe2 => probe3, but not probe1 => probe2 => probe3
            finishm_graph.probe_nodes.each_with_index do |node, probe_index2|
              # Don't add the node itself. This is a special case which is already handled below
              unless probe_index == probe_index2 or node.nil?
                terminal_nodes << finishm_graph.velvet_oriented_node(probe_index2).to_settable
              end
            end

            stack = DS::Stack.new
            stack.push finishm_graph.initial_path_from_probe(probe_index)

            while current_path = stack.pop
              log.debug "Perhaps #{current_path}?" if log.debug?
              if golden_onodes.include?(current_path.last.to_settable)
                # Probably a golden fragment, unless the node found is in the current path.
                # if that is true, that's a loop along a golden path.
                first_index = nil
                current_path.each_with_index do |directed_node, i|
                  if directed_node.node == current_path.last.node and
                    directed_node.first_side == current_path.last.first_side
                    first_index = i
                    break
                  end
                end

                if first_index == current_path.length-1
                  # Found another golden path(s)
                  log.debug "Ran into a golden node" if log.debug?
                  golden_path_fragments.push current_path
                  current_path.each do |onode|
                    golden_onodes << onode.to_settable
                  end
                else
                  # Loop found along a golden path or fragment
                  log.debug "Found a loop along a golden path: #{current_path}" if log.debug?
                  next
                end
              elsif already_visited_nodes.include?(current_path.last.to_settable) and
                current_path.last.node.node_id != terminal_node.node_id
                # Already seen this node, do nothing with it
                log.debug "Skipping #{current_path.last} since that has already been seen" if log.debug?
                next
              else
                if log.debug?
                  second_last_node = current_path[current_path.length-2]
                  second_last_node ||= 'initial_node'
                  log.debug "That's a new node, #{second_last_node}/#{current_path.last}" if log.debug?
                end

                # if we aren't beyond the leash. Presumably this
                # doesn't happen much, but there is a possibility that the leash will
                # prevent a real path being discovered. If there is two or more paths to a node
                # and a path longer than the leash is discovered first, then all the
                # nodes on that leash will be marked as discovered when they really aren't
                # TODO: fix the above bug
                if leash_length.nil? or current_path.length_in_bp < leash_length
                  # Found a new node for the user to play with
                  already_visited_nodes << current_path.last.to_settable

                  # Have we found a path to one of the other probes?
                  if terminal_nodes.include?(current_path.last.to_settable)
                    log.debug "Found the terminal node: #{current_path}" if log.debug?
                    golden_paths.push current_path

                  else # Use an else statement here because we want exploration to stop when other probes are encountered

                    # prep for next time
                    # Sort node IDs to simplify testing
                    next_nodes = current_path.neighbours_of_last_node(graph).sort{|n1, n2|
                      -(n1.node.node_id <=> n2.node.node_id)
                    }
                    next_nodes.each do |n|
                      path = current_path.copy
                      path.add_oriented_node n
                      log.debug "Pushing a path yet to be explored: #{path}" if log.debug?
                      stack.push path
                    end

                    # If we are at a dead end, add this whole path as a golden path. This
                    # is low coverage fluff, perhaps. Or it is nothing.
                    if next_nodes.empty?
                      log.debug "Found a dead end path: #{current_path}" if log.debug?
                      golden_paths.push current_path
                    end
                  end
                else
                  if log.debug?
                    log.debug "Found a path that made it to the end of the leash, with path length #{current_path.length_in_bp} vs leash length #{leash_length}"
                    log.debug "Path given up on: #{current_path}"
                    log.debug "Path sequence given up on: #{current_path.sequence}"
                    log.debug "Node lengths: #{current_path.collect{|n| n.node.length_alone}.join(',')}"
                  end
                  # Record this past-leash-length path
                  golden_paths.push current_path
                end
              end
            end

            log.debug "Found #{golden_paths.length} golden paths and #{golden_fragments.length} golden fragments" if log.debug?
            fluff_half_result = FlufferHalfResult.new
            fluff_half_result.golden_paths = golden_paths
            fluff_half_result.golden_path_fragments = golden_fragments

            half_results.push fluff_half_result
          end
        end

        return half_results
      end

      def fluff_part2(half_results)
        all_all_paths = []

        half_results.each do |segment_half_result|
          # OK, so we've transformed the data into a state where there is
          # at least one path through the data
          # and tentacles hanging off various golden nodes.
          # Now separate out the paths and return the array.
          # First transform the data so it can be referenced by the end node
          terminal_golden_nodes_to_paths = {}
          segment_half_result.golden_path_fragments.each do |fragment|
            l = fragment.last.to_settable
            terminal_golden_nodes_to_paths[l] ||= []
            terminal_golden_nodes_to_paths[l].push fragment
          end
          # Next backtrack through the paths
          # Each path starts at the beginning and ends at a
          # golden node
          all_paths = []
          stack = DS::Stack.new
          # Push the golden path and all paths that finish at the last node
          segment_half_result.golden_paths.each do |golden_path|
            stack.push [golden_path, 0]
          end

          while array = stack.pop
            current_path = array[0]
            num_to_ignore = array[1]

            log.debug "Defragging #{current_path.to_s}, ignoring the last #{num_to_ignore} nodes" if log.debug?
            all_paths.push current_path

            # Iterate down this path, spawning new paths if there
            # are paths that intersect
            passed_nodes = []
            current_path.trail.reverse.each_with_index do |onode, i|
              unless i < num_to_ignore #ignore the last one(s) because they've already been handled
                frags = terminal_golden_nodes_to_paths[onode.to_settable]
                log.debug "Offshoots from #{onode}: #{frags.nil? ? '[]' : frags.collect{|f| f.collect{|n| n.node_id}.join(',')}.join(' and ')}" if log.debug?
                if frags
                  frags.each do |fragment|
                    log.debug "Using an offshoot: #{fragment.to_s}" if log.debug?
                    # The fragment extends from the beginning to the golden node,
                    # the current node. So create a new complete path,
                    # And push it to the stack.
                    new_golden = fragment.copy
                    log.debug "Adding #{new_golden.to_s} and #{passed_nodes.collect{|n| n.node_id}}" if log.debug?
                    passed_nodes.reverse.each_with_index do |onode, i|
                      new_golden.add_oriented_node onode
                    end
                    log.debug "Enqueueing: #{new_golden.to_s} ignoring the last #{i+1} nodes" if log.debug?
                    stack.push [new_golden, i+1]
                  end
                end
              end
              passed_nodes.push onode
            end
          end

          # All the paths are in an array, just a linear series of distinct paths
          all_all_paths.push all_paths
        end

        return all_all_paths
      end
    end
  end
end

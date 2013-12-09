require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AcyclicConnectionFinder
      include Bio::FinishM::Logging

      def find_trails_between_nodes(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        #find_trails_between_nodes_depth_first_search(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)

        initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        way = start_looking_off_the_end_of_the_first_node ?
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
        initial_path.add_node initial_node, way

        return find_all_trails_squid_way(graph, initial_path, terminal_node, leash_length)
      end

      def find_all_trails_squid_way(graph, path, terminal_node, leash_length)
        half_result = find_all_trails_squid_way_part1(graph, path, terminal_node, leash_length)

        log.debug "Golden path: #{half_result.golden_path.to_s}" if log.debug?
        log.debug "found #{half_result.golden_fragments.length} golden fragments: #{half_result.golden_fragments.collect{|g| g.to_s}.join("\n")}" if log.debug?
        return [] if half_result.golden_path.nil?

        if half_result.golden_fragments.length > 10
          log.error "Too many paths found, not enumerating them because this would take too long"
          return []
        end

        return find_all_trails_squid_way_part2(half_result)
      end

      def find_all_trails_squid_way_part1(graph, path, terminal_node, leash_length)
        # do a depth first search, but keep track of which nodes are known to
        # to make it to the end. Then at the end work backwards to collect the paths
        # and separate them all.

        golden_onodes = Set.new #Nodes known to lead to the final destination
        golden_path = nil #This is the first full path to be found
        golden_path_fragments = [] #paths to join up to the end
        already_visited_nodes = Set.new #Nodes that have already been explored

        stack = DS::Stack.new
        stack.push path.copy

        # Depth first searching
        # While there is paths that haven't been fully explored
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
              log.debug "That's a new node, #{second_last_node}/#{current_path.last}"
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

              # Have we found a path?
              if current_path.last.node.node_id == terminal_node.node_id
                log.debug "Found the terminal node: #{current_path}" if log.debug?
                if golden_path.nil?
                  # This is the first path to be found
                  golden_path = current_path.copy
                  current_path.each do |onode|
                    golden_onodes << onode.to_settable
                  end
                else
                  # Found a new path that only converges on the terminal node
                  golden_path_fragments.push current_path.copy
                  current_path.each do |onode|
                    golden_onodes << onode.to_settable
                  end
                end
              end

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
            else
              log.warn "Giving up on this path because this is beyond the leash, with length #{path.length_in_bp} vs leash length #{leash_length}"
              if log.debug?
                log.debug "Path given up on: #{path}"
                log.debug "Path sequence given up on: #{path.sequence}"
                log.debug "Node lengths: #{path.collect{|n| n.node.length_alone}.join(',')}"
              end
            end
          end
        end
        half_result = HalfSquidResult.new
        half_result.golden_path = golden_path
        half_result.golden_fragments = golden_path_fragments
        return half_result
      end

      def find_all_trails_squid_way_part2(half_result)
        golden_path = half_result.golden_path
        golden_path_fragments = half_result.golden_fragments

        # OK, so we've transformed the data into a state where there is
        # at least one path through the data
        # and tentacles hanging off various golden nodes.
        # Now separate out the paths and return the array.
        # First transform the data so it can be referenced by the end node
        terminal_golden_nodes_to_paths = {}
        golden_path_fragments.each do |fragment|
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
        stack.push [golden_path, 0]

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
        return all_paths
      end

      class HalfSquidResult
        attr_accessor :golden_path, :golden_fragments
      end
    end
  end
end

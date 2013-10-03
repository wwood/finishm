require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AcyclicConnectionFinder
      def log
        Bio::Log::LoggerPlus['finishm']
      end

      def find_trails_between_nodes(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        #find_trails_between_nodes_depth_first_search(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)

        initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        way = start_looking_off_the_end_of_the_first_node ?
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
        initial_path.add_node initial_node, way

        return find_all_trails_squid_way(graph, initial_path, terminal_node, leash_length)
      end

      def find_all_trails_between_nodes(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        way = start_looking_off_the_end_of_the_first_node ?
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
        initial_path.add_node initial_node, way

        return depth_first_search_recursive(graph, initial_path, terminal_node, leash_length)
      end

      def find_all_trails_squid_way(graph, path, terminal_node, leash_length)
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
          #log.debug "Perhaps #{current_path}?" if log.debug?
          if golden_onodes.include?(current_path.last.to_settable)
            # Found another golden path(s)
            log.debug "Ran into a golden node" if log.debug?
            golden_path_fragments.push current_path
            current_path.each do |onode|
              golden_onodes << onode.to_settable
            end
          elsif already_visited_nodes.include?(current_path.last.to_settable) and
            current_path.last.node.node_id != terminal_node.node_id
            # Already seen this node, do nothing with it
            log.debug "Skipping #{current_path.last} since that has already been seen" if log.debug?
            next
          else
            log.debug "That's a new node, #{current_path.last}" if log.debug?

            # if we aren't beyond the leash. Presumably this
            # doesn't happen much, but there is a possibility that the leash will
            # prevent a real path being discovered. If there is two or more paths to a node
            # and a path longer than the leash is discovered first, then all the
            # nodes on that leash will be marked as discovered when they really aren't
            # TODO: fix the above bug
            if leash_length.nil? or path.length_in_bp < leash_length
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
                stack.push path
              end
            else
              log.debug "Giving up because this is beyond the leash" if log.debug?
            end
          end
        end
        log.debug "Golden path: #{golden_path.to_s}" if log.debug?
        log.debug "found #{golden_path_fragments.length} golden fragments: #{golden_path_fragments.collect{|g| g.to_s}.join("\n")}" if log.debug?
        log.debug '                 ' if log.debug?
        return [] if golden_path.nil?

        # OK, so we've transformed the data into a state where there is
        # at least one path through the data (or it is nil if there is no path)
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

      # A recursive depth first search algorithm for exploring the graph. Returns an array
      # of paths between the current path and the terminal node
      def depth_first_search_recursive(graph, path, terminal_node, leash_length)
        paths_found = []
        log.debug "Depth-first searching from #{path.to_s}"

        neighbours = path.neighbours_of_last_node(graph)

        # First, explore the neighbours - are we there yet?
        neighbours.each do |neighbour|
          log.debug "Considering neighbour in loop 1: #{neighbour}"
          #TODO: be cogent of orientation of the terminal node. Analogous fix for  second loop below
          if neighbour.node == terminal_node
            #yey, found a path. Record this.
            log.info "Found a trail"
            path.add_oriented_node neighbour
            paths_found.push path.copy
            path.remove_last_node
          end
        end

        # Next do further exploration
        neighbours.each do |neighbour|
          log.debug "Considering neighbour in loop 2: #{neighbour}"
          # Ignore cyclical paths
          next if path.include_oriented_node?(neighbour)

          # Ignore paths that are finished
          next if neighbour.node == terminal_node

          # DFS the other nodes
          path.add_oriented_node neighbour
          if path.length_in_bp > leash_length
            log.debug "Exploration truncated as leash length has been exceeded, along: #{path.to_s}"
          else
            paths_found.push depth_first_search_recursive(graph, path, terminal_node, leash_length)
          end
          path.remove_last_node
        end

        return paths_found.flatten
      end

      def find_trails_between_nodes_depth_first_search(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        paths_found = []
        discovered_list = Set.new
        explored_list = Set.new
        known_edges = Set.new

        stack = DS::Stack.new

        path = Bio::Velvet::Graph::OrientedNodeTrail.new
        way = start_looking_off_the_end_of_the_first_node ?
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
        path.add_node initial_node, way
        stack.push path

        while current_path = stack.pop
          current_node = current_path.last.node
          log.debug "Just popped #{current_node.node_id}"
          if current_node == terminal_node
            log.info "Found a path between initial and terminal nodes!"
            paths_found.push current_path
          end

          # Visit all the adjacent nodes
          str = current_path.collect{|n| "#{n.node.node_id}_#{n.node.coverage.round}"}.join(' ')
          log.debug "Finding next neighbours of this trail: #{str}"
          finished_exploring = true
          neighbours = current_path.neighbours_of_last_node(graph)
          neighbours.each do |neighbour|
            edges = graph.get_arcs_by_node current_node, neighbour.node
            raise "dragons" if edges.length != 1
            edge = edges[0]
            log.debug "Considering neighbour #{neighbour.node.node_id}" if log.debug?

            if known_edges.include?(edge)
              log.debug "Already seen this edge, ignoring: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
              next
            end
            known_edges << edge

            discovered = discovered_list.include?(neighbour.node)
            explored = explored_list.include?(neighbour.node)
            if !discovered and !explored
              log.debug "Found a new edge to discover/explore: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
              discovered_list << neighbour
              new_path = current_path.copy
              new_path.add_node neighbour.node, neighbour.first_side
              log.debug "Adding new path to the stack: #{new_path.to_s}"
              stack.push new_path
              log.debug "Stack is now #{stack.size} in length"
              finished_exploring = false
            end
          end

          if finished_exploring
            log.debug "Finished exploring #{current_path.collect{|ori_node|ori_node.node.node_id}.join(',')}"
            explored_list << current_path.last.node
            #popped = stack.pop
          end
        end
        log.info "Found #{paths_found.length} paths, after exploring #{explored_list.length} nodes and #{known_edges.length} edges"
        return paths_found
      end

      # Find a single path through the graph, ignoring the fact that the graph is directed.
      # Only useful for debug purposes, I imagine
      def find_trails_between_nodes_depth_first_search_undirected(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        paths_found = []
        discovered_list = Set.new
        explored_list = Set.new
        known_edges = Set.new

        stack = DS::Stack.new

        path = Bio::Velvet::Graph::OrientedNodeTrail.new
        way = start_looking_off_the_end_of_the_first_node ?
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
        path.add_node initial_node, way
        stack.push path

        while current_path = stack.pop
          current_node = current_path.last.node
          log.debug "Just popped #{current_node.node_id}"
          if current_node == terminal_node
            log.info "Found a path between initial and terminal nodes!"
            paths_found.push current_path
          end

          # Visit all the adjacent nodes
          finished_exploring = true
          last = current_path.last.node
          neighbours = [
            graph.neighbours_off_end(last),
            graph.neighbours_into_start(last)
          ].flatten.uniq.collect do |n|
            o = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
            o.node = n
            o.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
            o
          end
          neighbours.each do |neighbour|
            edges = graph.get_arcs_by_node current_node, neighbour.node
            raise "dragons" if edges.length != 1
            edge = edges[0]
            log.debug "Considering neighbour #{neighbour.node.node_id}" if log.debug?

            if known_edges.include?(edge)
              log.debug "Already seen this edge, ignoring: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
              next
            end
            known_edges << edge

            discovered = discovered_list.include?(neighbour.node)
            explored = explored_list.include?(neighbour.node)
            if !discovered and !explored
              log.debug "Found a new edge to discover/explore: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
              discovered_list << neighbour
              new_path = current_path.copy
              new_path.add_node neighbour.node, neighbour.first_side
              log.debug "Adding new path to the stack: #{new_path.to_s}"
              stack.push new_path
              log.debug "Stack is now #{stack.size} in length"
              finished_exploring = false
            end
          end

          if finished_exploring
            log.debug "Finished exploring #{current_path.collect{|ori_node|ori_node.node.node_id}.join(',')}"
            explored_list << current_path.last.node
            #popped = stack.pop
          end
        end
        log.info "Found #{paths_found.length} paths, after exploring #{explored_list.length} nodes and #{known_edges.length} edges"
        return paths_found
      end

      # Perform a search of the graph starting at the initial node, and try to
      # get all the paths between the initial and terminal nodes. So that things don't get out
      # of control, two bounds on this problem are given:
      #
      # 1. No cyclical trails are followed
      # 2. The total trail length (in base pairs) is less than the leash_length
      #
      # Return an array of Trail objects
      def find_trails_between_nodes_super_dumb(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        successful_trails = []
        search_underneath_nodes = lambda do |trail, new_node|
          trail.add_node new_node

          if new_node == terminal_node
            # We have found a successful trail. Done.
            log.info "Successful trail found: #{trail.collect{|node| node.node_id}.join(',')}" if log and log.debug?
            successful_trails.push trail

          else
            # OK so there is more to go (or we have strayed off the golden paths)

            # Need all the neighbours of the current node.
            neighbours = nil
            if trail.last_added_node_dangling_side == :end
              neighbours = graph.neighbours_off_end(new_node)
            elsif trail.last_added_node_dangling_side == :start
              neighbours = graph.neighbours_into_start(new_node)
            else
              raise "Programming error when creating trails between two nodes"
            end

            if neighbours.empty?
              # Dead end path
              log.debug "Dead end path found: #{trail}" if log and log.debug?
            else
              # There is one or more paths from here
              neighbours.each_with_index do |next_node, i|
                if trail.include?(next_node)
                  log.warn "Cyclical trail found in graph, assuming that cyclical path is a dead end." if log and log.debug?
                elsif trail.nucleotide_length + next_node.length_alone > leash_length
                  # We've strayed an unacceptable length away from the thing, forget it
                  log.debug "Straying too far trying to find paths, giving up on this trail. Length with next node is #{trail.nucleotide_length + next_node.length_alone}" if log and log.debug?
                else
                  # Keep on truckin'
                  log.debug "Truckin' onto the next node, from #{new_node.node_id} to #{next_node.node_id}" if log and log.debug?
                  arcs = graph.get_arcs_by_node new_node, next_node
                  arcs.reject!{|arc| arc.begin_node_id == arc.end_node_id}
                  if arcs.length != 1
                    raise "Programming error"
                  end
                  arc = arcs[0]
                  if arc.connects_to_beginning?(next_node.node_id)
                    trail.last_added_node_dangling_side = :end
                  elsif arc.connects_to_end?(next_node.node_id)
                    trail.last_added_node_dangling_side = :start
                  else
                    raise "Programming error"
                  end
                  search_underneath_nodes.call(trail.copy, next_node)
                end
              end
            end
          end
          log.debug "Ending recursion" if log and log.debug?
        end

        # initialise the search
        current_trail = Trail.new
        if start_looking_off_the_end_of_the_first_node
          current_trail.last_added_node_dangling_side = :end
        else
          current_trail.last_added_node_dangling_side = :start
        end
        search_underneath_nodes.call current_trail, initial_node

        return successful_trails
      end

      # Find trails between a set of different nodes
      def find_trails_between_node_set(graph, anchoring_nodes_and_directions, leash_length)
        path_sets = []

        # all vs all searching, but only the top triangle of the matrix
        anchoring_nodes_and_directions.each_with_index do |pair, i|
          initial_node = anchoring_nodes_and_directions.to_a[i][0]
          initial_node_direction = anchoring_nodes_and_directions.to_a[i][1]
          path = Bio::Velvet::Graph::OrientedNodeTrail.new
          way = initial_node_direction ?
            Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
            Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
          path.add_node initial_node, way

          target_nodes = Set.new
          j = 0
          anchoring_nodes_and_directions.each do |node, dir|
            #target all nodes that come after this one in the array
            target_nodes << node unless j<=i
            j += 1
          end

          path_sets.push find_trails_between_one_node_and_a_set_of_others(graph, path, target_nodes, leash_length)
        end
        return path_sets
      end

      def find_trails_between_one_node_and_a_set_of_others(graph, initial_path, target_nodes, leash_length)
        paths_found = []
        discovered_list = Set.new
        explored_list = Set.new
        known_edges = Set.new
        stack = DS::Stack.new

        path = initial_path.copy
        stack.push path

        while current_path = stack.pop
          current_node = current_path.last.node
          log.debug "Just popped #{current_node.node_id}"
          if target_nodes.include?(current_node)
            log.info "Found a path between two nodes nodes: #{initial_path.first.node.node_id} and #{current_node.node_id}!"
            paths_found.push current_path
          end

          # Visit all the adjacent nodes
          str = current_path.collect{|n| "#{n.node.node_id}_#{n.node.coverage.round}"}.join(' ')
          log.debug "Finding next neighbours of this trail: #{str}"
          finished_exploring = true
          neighbours = current_path.neighbours_of_last_node(graph)
          neighbours.each do |neighbour|
            edges = graph.get_arcs_by_node current_node, neighbour.node
            edges.each do |edge|
              log.debug "Considering neighbour #{neighbour.node.node_id}" if log.debug?

              if known_edges.include?(edge)
                log.debug "Already seen this edge, ignoring: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
                next
              end
              known_edges << edge

              discovered = discovered_list.include?(neighbour.node)
              explored = explored_list.include?(neighbour.node)
              if !discovered and !explored
                log.debug "Found a new edge to discover/explore: #{edge.begin_node_id}/#{edge.end_node_id}" if log.debug?
                discovered_list << neighbour
                new_path = current_path.copy
                new_path.add_node neighbour.node, neighbour.first_side
                if new_path.length_in_bp > leash_length
                  log.debug "Discontinuing this path because it is beyond the leash length: #{new_path.to_s}"
                else
                  log.debug "Adding new path to the stack: #{new_path.to_s}"
                  stack.push new_path
                  log.debug "Stack is now #{stack.size} in length"
                  finished_exploring = false
                end
              end
            end
          end

          if finished_exploring
            log.debug "Finished exploring #{current_path.collect{|ori_node|ori_node.node.node_id}.join(',')}"
            explored_list << current_path.last.node
            #popped = stack.pop
          end
        end
        log.info "Found #{paths_found.length} paths, after exploring #{explored_list.length} nodes and #{known_edges.length} edges"
        return paths_found
      end

      class Trail < Array
        attr_accessor :last_added_node_dangling_side

        attr_reader :length_in_bp

        def add_node(node)
          push node
          @length_in_bp ||= 0
          @length_in_bp += node.length_alone
        end

        # Number of nucleotides in this trail
        def nucleotide_length
          return 0 if empty?
          len = self[0].length
          each_with_index do |node, i|
            next if i == 0
            len += node.length_alone
          end
          return len
        end

        # Return a copy of the current trail such that modifying the path
        # of the new node does not affect the path of the original
        def copy
          new_trail = Trail.new
          new_trail.last_added_node_dangling_side = @last_added_node_dangling_side
          each do |node|
            new_trail.push node
          end
          return new_trail
        end

        def to_s
          "Trail #{self.object_id}: #{collect{|n| n.node_id}.join(',')}"
        end
      end
    end
  end
end

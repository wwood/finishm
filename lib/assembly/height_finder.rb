require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::HeightFinder
  include Bio::FinishM::Logging

  # visit nodes in range and determine heights
  def traverse(graph, options={})
    range = options[:range]
    forwards = !options[:reverse]
    by_height = []
    traversal_nodes = {}
    nodes_in_retrace_phase = Set.new
    if range.nil?
      nodes = []
      graph.nodes.each do |node|
        nodes.push Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new node, forwards
      end
    else
      nodes = range
    end

    # depth-first so stack
    stack = DS::Stack.new
    nodes.each do |onode|
      traversal_node = TraversalNode.new
      traversal_node.onode = onode
      traversal_node.nodes_in = []
      traversal_nodes[onode.to_settable] = traversal_node
      stack.push traversal_node
    end

    while traversal_node = stack.pop
      settable = traversal_node.onode.to_settable
      describe = nil

      if log.debug?
        describe = "node #{settable[0]} direction #{settable[1]}"
        log.debug "visiting #{describe}."
      end

      solved = !traversal_node.height.nil?
      if solved
        log.debug "Height of #{describe} is already known. Skip." if log.debug?
        next
      end

      # find neighbours
      neighbours = traversal_node.nodes_out
      if neighbours.nil?
        neighbours = traversal_node.onode.next_neighbours(graph)
        unless range.nil?
          neighbours.reject!{|onode| range.include? onode} #not in defined range
        end

        # Get or create traversal version of node
        neighbours = neighbours.collect do |onode|
          nbr_settable = onode.to_settable
          traversal_nbr = traversal_nodes[nbr_settable]
          if traversal_nbr.nil?
            traversal_nbr = TraversalNode.new
            traversal_nbr.onode = onode
            traversal_nbr.nodes_in = []
            traversal_neighbours[nbr_settable] = traversal_nbr
          end
          traversal_nbr
        end

        #remember neighbours
        traversal_node.nodes_out = neighbours
      end


      # Can we solve the node?
      if neighbours.empty? #check for a tip
        log.debug "#{describe} is a tip." if log.debug?
        traversal_node.height = 0
        if by_height[0].nil?
          by_height[0] = [traversal_node]
        else
          by_height[0].push(traversal_node)
        end
        log.debug "Found height of #{describe}: 0." if log.debug?
        next
      end

      if nodes_in_retrace_phase.include? settable
        log.debug "Retracing back to #{describe}." if log.debug?
        # Neighbours should have been explored, unsolved neighbours implies a closed cyclic path
        solved_neighbours = neighbours.reject{|node| node.height.nil?}
        unless solved_neighbours.empty?
          log.debug "We know the heights of neighbours #{solved_neighbours.collect{|node| node.onode.node.node_id}.join(',')}." if log.debug?
          # Compute height from solved neighbours
          height = solved_neighbours.map{|node| node.height}.max + 1
          log.debug "Found height of #{describe}: #{height}." if log.debug?
          traversal_node.height = height
          if by_height[height].nil?
            by_height[height] = [traversal_node]
          else
            by_height[height].push(traversal_node)
          end
        end
        # If no solved neighbours, leave unsolved (node joins closed cyclic path)

        # Move out of retrace phase
        nodes_in_retrace_phase.delete settable
        log.debug "Finished retracing #{describe}." if log.debug?
        next
      end

      # Return node stack marked explored and push neighbours
      nodes_in_retrace_phase << settable
      stack.push traversal_node
      log.debug "Pushing #{describe} in retrace mode." if log.debug?
      neighbours.each do |node|
        node_settable = node.onode.to_settable

        # Note the parent of neighbour unless already known
        nodes_in = node.nodes_in
        if not nodes_in.any?{|nbr| nbr.onode.to_settable == node_settable}
          nodes_in.push traversal_node
        end

        if nodes_in_retrace_phase.include? node_settable
          # A currently retracing neighbour implies a cycle, cut it off here
          log.debug "Neighbour node #{node_settable[0]} direction #{node_settable[1]} is retracing. Not revisiting." if log.debug?
        else
          log.debug "Pushing neighbour node #{node_settable[0]} direction #{node_settable[1]}." if log.debug?
          stack.push node
        end
      end
    end
    return by_height
  end

  # maximum paths
  def max_paths_through(by_height)
    max_paths_from = {}
    by_height.each_with_index do |nodes, level|
      log.debug "At height #{level}." if log.debug?
      if level == 0  # tips
        nodes.each do |node|
          settable = node.onode.to_settable
          log.debug "Counted maximum of 1 path to node #{settable[0]} direction #{settable[1]}." if log.debug?
          max_paths_from[settable] = 1
        end
        next
      end

      nodes.each do |node|
        settable = node.onode.to_settable
        max_paths_from_neighbours = node.nodes_out.collect{|nbr| max_paths_from[nbr.onode.to_settable]}.reject{|n| n.nil?}
        log.debug "Found neighbours of node #{settable[0]} direction #{settable[1]} with maximum paths #{max_paths_from_neighbours.join(',')}." if log.debug?
        max_paths_from[settable] = max_paths_from_neighbours.reduce{|memo, num| memo+num}
        log.debug "Counted maximum of #{max_paths_from[settable]} paths to node #{settable[0]} direction #{settable[1]}." if log.debug?
      end
    end

    # Get the graph roots (which are nodes with no parents) and add max_paths_from for each to get graph total
    if log.debug?
      all_nodes = by_height.flatten
      log.debug "all nodes #{all_nodes.collect{|node| node.onode.node_id}.join(',')}."
      nodes_without_parents = all_nodes.select{|node| node.nodes_in.empty?}
      log.debug "nodes with no parents #{nodes_without_parents.collect{|node| node.onode.node_id}.join('.')}."
      root_settables = nodes_without_parents.collect{|node| node.onode.to_settable}
      log.debug "settables for roots #{root_settables}"
    end
    root_keys = by_height.flatten.select{|node| node.nodes_in.empty? }.collect{|node| node.onode.to_settable}
    log.debug "Found graph roots #{root_keys.collect{|settable| settable[0]}.join(',')} with maximum paths #{root_keys.collect{|key| max_paths_from[key]}.join(',')}." if log.debug?
    max_paths = root_keys.map{|settable| max_paths_from[settable]}.reduce{|memo, num| memo+num}
    log.debug "Counted maximum of #{max_paths} through graph." if log.debug?
    return max_paths
  end

  # minimum paths
  def min_paths_through(by_height)
    live_nodes = Set.new
    max_alive_counter = 0
    by_height.each_with_index do |nodes, level|
      log.debug "At height #{level}." if log.debug?
      # nodes at current level become live
      nodes.each do |node|
        settable = node.onode.to_settable
        log.debug "Setting node #{settable[0]} direction #{settable[1]} as live." if log.debug?
        live_nodes << settable
      end
      if level > 0
        #children of nodes at current level are no longer live
        nodes.each do |node|
          children = node.nodes_out
          children.each do |node|
            settable = node.onode.to_settable
            log.debug "Setting child node #{settable[0]} direction #{settable[1]} of live node #{node.onode.node_id} as inactive." if debug?
            live_nodes.delete(node.onode.to_settable)
          end
        end
      end

      log.debug "There are currently #{live_nodes.length} nodes alive. Max is #{max_alive_counter}." if log.debug?
      if live_nodes.length > max_alive_counter
        #track the maximum live nodes at any level
        log.debug "Updating max to #{live_nodes.length}." if log.debug?
        max_alive_counter = live_nodes.length
      end
    end
    return max_alive_counter
  end

  class TraversalNode
    attr_accessor :onode, :height, :nodes_in, :nodes_out
  end

  # bubble finder
  class BubbleFinder
    def find_bubbles(height_finder)
      num_live_forks_from_node = {}
      seen_forks = {}
      bubbles = {}
      height_finder.ascend_by_level do |nodes, level|
        if level == 0
          # tips
          nodes.each do |onode|
            settable = onode.to_settable
            if height_finder.nodes_in.has_key? settable and height_finder.nodes_in[settable].length > 1
              num_live_forks_from_node[settable] = height_finder.nodes_in[settable].length
              fork_trail_set = ForkTrailSet.new
              fork_trail = ForkTrailSet::ForkTrail.new
              fork_trail.onodes = [onode]
              fork_trail_set.add fork_trail
              seen_forks[settable] = fork_trail_set
              bubbles[settable] = Set.new
              bubbles[settable] << settable
            end
          end
          next
        end

        nodes.each do |onode|
          settable = onode.to_settable
          child_forks = height_finder.nodes_out[settable].collect{|onode| seen_forks[onode.to_settable]}.reject{|f| f.nil?}
          fork_trail_set = nil
          if not child_forks.empty?
            # merge all child fork sets
            fork_trail_set = ForkTrailSet.new
            child_forks.each{|forks| fork_trail_set.merge forks}

            # register bubble membership
            member_of_bubble = Set.new
            fork_trail_set.all_forks{|onode| memeber_of_bubble << onode.to_settable}
            member_of_bubble.each{|bubble_settable| bubbles[bubble_settable] << settable}

            # deduplicate and try to close bubbles
            to_join = fork_trail_set.deduplicate

            while trail_to_join = to_join.pop
              onode = trail_to_join.onodes.last
              num_live_forks_from_node[onode.to_settable] -= 1
              if num_live_forks_from_node[onode.to_settable] == 0
                fork_trail_set.close trail_to_join
                fork_trail_set.deduplicate.each{|trail| to_join.push trail}
              end
            end
          end

          # fork away!
          if height_finder.nodes_in.has_key? settable and height_finder.nodes_in[settable].length > 1
            num_live_forks_from_node[settable] = height_finder.nodes_in[settable].length
            if fork_trail_set.nil?
              fork_trail_set = ForkTrailSet.new
              fork_trail = ForkTrailSet::ForkTrail.new
              fork_trail.onodes = [onode]
              fork_trail_set.push fork_trail
            else
              fork_trail_set.push onode
            end
          end
          seen_forks[settable] = fork_trail_set
        end
      end
      return bubbles.values
    end

    class ForkTrailSet
      attr_accessor :fork_trails

      def initialize
        @fork_trails = []
      end

      def all_forks
        @fork_trails.each{|trail| trail.onodes.each{|onode| yield onode}}
      end

      def deduplicate
        seen_nodes = Set.new
        duplicates = []
        keep = []
        @fork_trails.each do |trail|
          onode = trail.onodes.last
          if seen_nodes.include? onode.to_settable
            duplicates.push trail
          else
            keep.push trail
            seen_nodes << onode.to_settable
          end
        end
        @fork_trails = keep
        return duplicates
      end

      def merge(fork_trail_set)
        fork_trail_set.fork_trails.each{|trail| @fork_trails.push trail.copy}
      end

      def push(onode)
        @fork_trails.each{|trail| trail.onodes.push onode}
      end

      def close(fork_trail)
        @fork_trails.each do |trail|
          if trail.onodes.last == fork_trail.onodes.last
            trail.onodes.pop
          end
        end
        @fork_trails.reject!{|trail| trail.onodes.empty?}
      end

      class ForkTrail
        attr_accessor :onodes

        def matches?(fork_trail)
          return @onodes.last == fork_trail.onodes.last
        end

        def copy
          fork_trail = ForkTrail.new
          fork_trail.onodes = @onodes[0..-1]
        end
      end
    end
  end
end

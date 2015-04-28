require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::HeightFinder
  include Bio::FinishM::Logging

  attr_reader :nodes_in, :nodes_out, :by_height

  def initialize(graph, range = nil)
    @graph = graph
    @range = range
  end

  # visit nodes in range and determine heights
  def traverse()
    @by_height = []
    @nodes_in = {}
    @nodes_out = {}
    node_height = {}
    nodes_in_retrace_phase = Set.new
    nodes = @range || @graph.nodes

    # depth-first so stack
    stack = DS::Stack.new
    nodes.each{|onode| stack.push onode}

    while node = stack.pop
      settable = node.to_settable
      describe = nil

      if log.debug
        describe = "node ${settable[0]} direction ${settable[1]}."
        log.debug "visiting ${describe}."
      end

      solved = node_height.has_key? settable
      if solved
        log.debug "Height of ${describe} is already known. Skip." if log.debug?
        next
      end

      # find neighbours
      if @nodes_out.has_key? settable
        neighbours = @nodes_out[settable]
      else
        neighbours = node.next_neighbours(@graph)
        unless range.nil?
          neighbours.reject!{|onode| range.include? onode} #not in defined range
        end
        @nodes_out[settable] = neighbours #remember neighbours
      end


      # Can we solve the node?
      if neighbours.empty? #check for a tip
        log.debug "${describe} is a tip." if log.debug?
        node_height[settable] = 0
        @by_height[0] |= []
        @by_height[0].push(node)
        next
      end

      if nodes_in_retrace_phase.include? settable
        # Neighbours should have been explored, unsolved neighbours implies a closed cyclic path
        solved_neighbours = neighbours.select{|onode| node_height.has_key? onode }
        unless solved_neighbours.empty?
          # Compute height from solved neighbours
          height = solved_neighbours.map{|onode| node_height[onode.to_settable]}.max + 1
          node_height[settable] = height
          @by_height[height] |= []
          @by_height[height].push(node)
        end
        # If no solved neighbours, leave unsolved (node joins closed cyclic path)

        # Move out of retrace phase
        nodes_in_retrace_phase.delete settable
        next
      end

      # Return node stack marked explored and push neighbours
      nodes_in_retrace_phase << settable
      stack.push node
      neighbours.each do |onode|
        @nodes_in[onode.to_settable] |= Set.new
        @nodes_in[onode.to_settable] << settable
        unless nodes_in_retrace_phase.include? onode.to_settable
          # A currently retracing neighbour implies a cycle, cut it off here
          stack.push onode
        end
      end
    end
    @traversed = true
  end

  def ascend_by_level
    raise if not @traversed
    # Ascend heights to compute paths
    @by_height.each_with_index
  end

  # maximum paths
  def max_paths_through
    max_paths_from = {}
    ascend_by_level do |nodes, level|
      if level == 0  # tips
        nodes.each{|onode| max_paths_from[onode.to_settable] = 1}
        next
      end

      nodes.each do |onode|
        settable = onode.to_settable
        max_paths_from_neighbours = @nodes_out[settable].collect{|onode| max_paths_from[settable]}
        max_paths_from[settable] = max_paths_from_neighbours.reduce{|memo, num| memo+num}
      end
    end

    # Get the graph roots (which are nodes with no parents) and add max_paths_from for each to get graph total
    root_keys = Set.new( (@range || @graph.nodes).collect{|onode| onode.to_settable}.reject{|settable| @nodes_in.include? settable} )
    return max_paths_from.select{|key, val| root_keys.include? key}.values.reduce{|memo, num| memo+num}
  end

  # minimum paths
  def min_paths_through
    live_nodes = Set.new
    max_alive_counter = 0
    ascend_by_level do |nodes, level|
      # nodes at current level become live
      nodes.each{|onode| live_nodes << onode.to_settable}
      if level > 0
        #children of nodes at current level are no longer live
        nodes.each do |onode|
          children = @nodes_out[onode.to_settable]
          children.each{|onode| live_nodes.delete(onode.to_settable)}
        end
      end
      if live_nodes.length > max_path_counter
        #track the maximum live nodes at any level
        max_alive_counter = live_nodes.length
      end
    end
    return max_alive_counter
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

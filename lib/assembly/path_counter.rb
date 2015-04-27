require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::PathCounter
  include Bio::FinishM::Logging

  def initialize(graph, range = nil)
    @graph = graph
    @range = range
  end

  #
  def traverse()
    @by_height = []
    @nodes_in = {}
    @nodes_out = {}
    node_height = {}
    explored_nodes = Set.new
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

      explored = explored_nodes.include? settable
      solved = node_height.has_key? settable

      if solved
        log.debug "Height of ${describe} is already known. Skip." if log.debug?
        next
      end

      # find neighbours
      neighbours = node.next_neighbours(@graph)
      unless range.nil?
        neighbours.reject!{|onode| range.include? onode} #not in defined range
      end
      @nodes_out[settable] = neighbours #remember neighbours


      # Can we solve the node?
      if neighbours.empty? #check for a tip
        log.debug "${describe} is a tip." if log.debug?
        node_height[settable] = 0
        @by_height[0] |= []
        @by_height[0].push(node)
        explored_nodes << settable
        next
      end

      if explored
        # neighbours should be explored, unsolved neighbours are probably involved in cycles
        unless neighbours.collect{|onode| explored_nodes.include? onode.to_settable }.all?
          log.debug "Unexplored neighbours of an explored node???" if log.debug?
          raise
        end

        solved_neighbours = neighbours.select{|onode| node_height.has_key? onode }
        unless solved_neighbours.empty?
          # Compute height from solved neighbours
          height = solved_neighbours.map{|onode| node_height[onode.to_settable]}.max + 1
          node_height[settable] = height
          @by_height[height] |= []
          @by_height[height].push(node)
        end

        #If no solved neighbours, leave unsolved (probably involved in a cycle)
        #Cycles are tricky :(
        next
      end

      #Return node stack marked explored and push neighbours
      explored_nodes << settable
      stack.push node
      neighbours.each do |onode|
        @nodes_in[onode.to_settable] |= []
        @nodes_in[onode.to_settable].push(node)
        stack.push onode
      end
    end
    @traversed = true
  end

  def backtrack(&block)
    # Ascend heights to compute paths
    @by_height.each{|nodes| nodes.each{|onode| block.call(onode)}}
  end

  # maximum paths
  def max_paths_from_nodes
    max_paths_from = {}
    backtrack do |onode|
      settable = onode.to_settable
      unless @nodes_out.has_key? settable
        max_paths_from[settable] = 1
        next
      end
      max_paths_from_neighbours = @nodes_out[settable].collect{|onode| max_paths_from[settable]}
      max_paths_from[settable] = max_paths_from_neighbours.reduce{|memo, num| memo+num}
    end
    return max_paths_from
  end

  # minimum paths
  def min_paths_from_nodes
    min_paths_from = {}
    seen_forks = {}
    backtrack do |onode|
      settable = onode.to_settable
      unless @nodes_out.has_key? settable
        if @nodes_in.has_key? settable and @nodes_in[settable].length > 1
          bubbling[settable] = bubble
        end
        min_paths_from[settable] = 1
        next
      end

      bubbles_from_neighbours = @nodes_out[settable].collect{|onode| bubbling[onode.to_settable]}.reject{|bubble| bubble.nil?}
      if not bubbles_from_neighbors.empty?
        # try to close bubbles
        bubble = BubbleFinder.new
        bubble.join bubbles_from_neighbours
      end

      if @nodes_in.has_key? settable and @nodes_in[settable].length > 1
        bubble = BubbleFinder.new
      end
      forks_from_neighbours = @nodes_out[settable].collect{|onode| @nodes_in[onode.to_settable].length}
    end
  end

  class ForkNode
    attr_accessor :count, :fork_trails

    def initialize
      @fork_trails = []
    end

    def matches(fork_trail)
      matches = []
      @fork_trails.each do |trail|
        if trail.matches? fork_trail
          matches.push trail.onodes.last
        end
      end
      return matches
    end

    def close(onode)
      @fork_trails.each do |trail|
        if trail.onodes.last == onode
          trail.pop
        end
      end
    end

    def join(fork_trail)
      @fork_trails.push fork_trail
    end
  end

  class ForkTrail
    include Enumerable

    attr_accessor :onodes

    def matches?(fork_trail)
      return onodes.last == fork_trail.onodes.last
    end

    def trim_match(fork_trail)
      trimmed = nil
      if matches? fork_trail
        trimmed = @onodes.pop
        fork_trail.onodes.pop
      end
      return trimmed
    end

    def copy
      fork_trail = ForkTrail.new
      fork_trail.onodes = @onodes[0..-1]
    end
  end
end

require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::Dijkstra
  include Bio::FinishM::Logging

  # Return an array of DistancedOrientedNode objects, those reachable from
  # the initial_oriented_node. options:
  # :leash_length => max distance explored,
  #    can be set to nil to search indefinitely
  # :ignore_directions: => true or false (default). If true, explore direction-independently.
  #   i.e. if 1s->3s and 2s->3s, then include 2s in the returned set of min_distances
  #   and continue exploring from 2s. Return each found node twice, once for each direction
  # :neighbour_finder => an object that responds to #neighbours(oriented_node) and
  #   returns an array of Bio::FinishM::PairedEndNeighbourFinder::Neighbour objects
  #   default: just search using OrientedNode#next_neighbours
  # :max_nodes => maximum number of nodes to return, to prevent out of control
  #   exploring of the graph. If there is plenty of nodes to explore, then the
  #   length of the returned hash is options[:max_nodes]+1 (+1 because the starting
  #   node is included). It will probably be longer if :ignore_directions == true, in that
  #   case the number of node_ids is constrained. It may also be longer if there is ties
  #   at the edges of the constrained exploration.
  #
  # Returns a Hash of [node_id, first_side] => distance
  def min_distances(graph, initial_oriented_node, options={})
    pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
    first = DistancedOrientedNode.new
    first.node = initial_oriented_node.node
    first.first_side = initial_oriented_node.first_side
    first.distance = 0
    pqueue.push first, first.distance

    to_return = {}
    first_node = true
    found_nodes = Set.new([first.node.node_id])

    while min_distanced_node = pqueue.shift

      # Add/overwrite the current one
      to_return[min_distanced_node.to_settable] = min_distanced_node.distance

      log.debug "Working from #{min_distanced_node.inspect}" if log.debug?

      if options[:leash_length] and min_distanced_node.distance > options[:leash_length]
        # we are passed leash length, and this is the nearest node. So we are finito.
        log.debug "passed the leash length, cutting short our travels" if log.debug?
        break
      end

      if options[:max_nodes] and found_nodes.length > options[:max_nodes]
        log.debug "passed max-nodes threshold and have #{found_nodes.length} nodes" if log.debug?
        # remove extras that may have been queued if we are over the limit
        distances_direction_agnostic = {}
        to_return.each do |key, distance|
          prev = distances_direction_agnostic[key[0]]
          if prev.nil? or prev > distance
            distances_direction_agnostic[key[0]] = distance
          end
        end
        if distances_direction_agnostic.length > options[:max_nodes]
          sorted = distances_direction_agnostic.to_a.sort{|a,b| a[1]<=>b[1]}
          # deal with ties i.e. at the edge there can be multiple neighbours
          last_distance = sorted[options[:max_nodes]][1]

          # only keep those nodes that are sufficiently close
          to_return.select! do |key, distance|
            distance <= last_distance
          end
        end
        break
      end

      # Queue nodes after this one
      current_distance = min_distanced_node.distance

      # Find neighbouring nodes
      neighbours = nil
      if options[:neighbour_finder]
        neighbours = options[:neighbour_finder].neighbours(min_distanced_node)
      else
        neighbours = min_distanced_node.next_neighbours(graph)
      end

      # explore each neighbour node
      neighbours.each do |onode|
        found_nodes << onode.node.node_id
        new_distance = current_distance
        if options[:neighbour_finder]
          # Don't use negative distances as this algorithm cannot handle it, and it is impossible
          # anyway
          if onode.distance > 0
            new_distance += onode.distance
          else
            new_distance += 0
          end
        end
        unless first_node
          new_distance += min_distanced_node.node.length_alone
        end

        if to_return[onode.to_settable] and to_return[onode.to_settable] <= new_distance
          # We already know a shorter path to this neighbour, so ignore it
          log.debug "Already seen this node at the same or shorter distance, going no further" if log.debug?
        else
          log.debug "Queuing new distance for neighbour: #{onode}: #{new_distance}" if log.debug?
          # new shortest distance found. queue it up
          distanced_node = DistancedOrientedNode.new
          distanced_node.node = onode.node
          distanced_node.first_side = onode.first_side
          distanced_node.distance = new_distance
          to_return[onode.to_settable] = new_distance
          pqueue.push distanced_node, new_distance

          if options[:ignore_directions]
            reverse = DistancedOrientedNode.new
            reverse.node = onode.node
            reverse.first_side = onode.reverse.first_side
            reverse.distance = new_distance
            to_return[onode.to_settable] = new_distance
            pqueue.push reverse, new_distance
          end
        end
      end

      first_node = false
    end

    # if ignore directions, then fixup the return so that each direction is included
    if options[:ignore_directions]
      new_to_return = {}
      to_return.each do |key, distance|
        keys = [
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST,
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST].collect do |direction|
            [key[0], direction]
          end
        new_distance = keys.collect{|k| to_return[k]}.reject{|d| d.nil?}.min
        keys.each do |key|
          new_to_return[key] = new_distance
        end
      end
      to_return = new_to_return
    end

    return to_return
  end

  # like #min_distances except explores in both directions
  def min_distances_in_both_directions(graph, node, options={})
    all_min_distances = {}
    [
      Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST,
      Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST,
      ].each do |direction|
        onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(node, direction)
        min_distances = min_distances(graph, onode, options)
        min_distances.each do |node_direction, distance|
          current = all_min_distances[node_direction]
          unless current and current > distance
            all_min_distances[node_direction] = distance
          end
        end
      end
    return all_min_distances
  end

  def min_distances_from_many_nodes_in_both_directions(graph, nodes, options={})
    all_min_distances = {}
    nodes.each do |node|
    [
      Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST,
      Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST,
      ].each do |direction|
        onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(node, direction)
        min_distances = min_distances(graph, onode, options)
        min_distances.each do |node_direction, distance|
          current = all_min_distances[node_direction]
          unless current and current > distance
            all_min_distances[node_direction] = distance
          end
        end
      end
    end
    return all_min_distances
  end

  # Generate an array of OrientedTrails
  def all_trails(graph, initial_oriented_node, options={})
    pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
    first = DistancedOrientedNodeTrail.new
    first_trail = OrientedNodeTrail.new
    first_trail.add_oriented_node(initial_oriented_node.copy)
    first.oriented_trail = first_trail
    first.distance = 0
    pqueue.push first, first.distance

    to_return = []
    first_node = true

    while min_distanced_trail = pqueue.shift

      # Add/overwrite the current one
      to_return.push min_distanced_trail

      log.debug "Working from #{min_distanced_trail.inspect}" if log.debug?

      if options[:leash_length] and min_distanced_trail.distance > options[:leash_length]
        # we are passed leash length, and this is the nearest node. So we are finito.
        # TIM - don't push more neighbours but finish current queue
        log.debug "passed the leash length, cutting short our travels" if log.debug?
        next
      end

      if options[:max_nodes] and found_nodes.length > options[:max_nodes]
        log.debug "passed max-nodes threshold and have #{found_nodes.length} nodes" if log.debug?
        # remove extras that may have been queued if we are over the limit
        next
      end

      # Queue nodes after this one
      current_distance = min_distanced_trail.distance

      # Find neighbouring nodes
      neighbours = nil
      if options[:neighbour_finder]
        neighbours = options[:neighbour_finder].neighbours(min_distanced_trail.last)
      else
        neighbours = min_distanced_trail.neighbours_of_last_node(graph)
      end

      # explore each neighbour node
      neighbours.each do |onode|
        found_nodes << onode.node.node_id
        new_distance = current_distance
        if options[:neighbour_finder]
          # Don't use negative distances as this algorithm cannot handle it, and it is impossible
          # anyway
          if onode.distance > 0
            new_distance += onode.distance
          else
            new_distance += 0 # TIM-???
          end
        end
        unless first_node
          new_distance += min_distanced_trail.last.node.length_alone
        end

        #Queue everything!!!
        log.debug "Queuing trail with distance to neighbour: #{onode}: #{new_distance}" if log.debug?
        distanced_trail = DistancedOrientedNodeTrail.new
        distanced_trail.trail = min_distanced_trail.copy
        distanced_trail.trail.add_oriented_node onode
        distanced_trail.distance = new_distance
        #to_return[onode.to_settable] = new_distance # TIM - Push to results when popped of queue
        pqueue.push distanced_trail, new_distance

        if options[:ignore_directions]
          reverse = DistancedOrientedTrail.new
          reverse.trail = min_distanced_trail.copy
          reverse.trail.trail.add_oriented_node onode.reverse
          reverse.distance = new_distance
          #to_return[onode.to_settable] = new_distance
          pqueue.push reverse, new_distance
        end
      end

      first_node = false
    end

    if options[:max_nodes] and found_nodes.length > options[:max_nodes]
      distances_direction_agnostic = {}
      to_return.each do |dtrail|
        last_node_id = dtrail.last.node.node_id
        prev = distances_direction_agnostic[last_node_id]
        if prev.nil? or prev > dtrail.distance
          distances_direction_agnostic[last_node_id] = dtrail.distance
        end
      end
      if distances_direction_agnostic.length > options[:max_nodes]
        sorted = distances_direction_agnostic.to_a.sort{|a,b| a[1]<=>b[1]}
        # deal with ties i.e. at the edge there can be multiple neighbours
        last_distance = sorted[options[:max_nodes]][1]

        # only keep those trails to nodes that are sufficiently short
        to_return.select! do |dtrail|
          dtrail.distance <= last_distance
        end
      end
    end

    # Return trails
    return to_return.map{|d| d.oriented_trail}
  end

  # Returns a Hash of [node_id, first_side] => [position_distance] => position_multiplicity
  def trail_stats(graph, initial_oriented_node, options={})
    counter = DistancedNodeCounter.new
    first = DistancedOrientedNode.new
    first.node = initial_oriented_node.node
    first.first_side = initial_oriented_node.first_side
    first.positions = Hash.new 1
    first.positions[0] = 1
    first.distance = 0
    counter.push first

    first_node = true
    found_nodes = Set.new([first.node.node_id])

    while min_distanced_node = counter.shift

      log.debug "Working from #{min_distanced_node.inspect}" if log.debug?
      log.debug "Node currently maps to the following distances: #{min_distanced_node.positions.to_a.collect{|d| d.join(': ')}.join(', ')}." if log.debug?

      if options[:leash_length] and min_distanced_node.distance > options[:leash_length]
        # we are passed leash length, and this is the nearest node. So we are finito.
        log.debug "passed the leash length, cutting short our travels" if log.debug?
        next
      end

      if options[:max_nodes] and found_nodes.length > options[:max_nodes]
        log.debug "passed max-nodes threshold and have #{found_nodes.length} nodes" if log.debug?
        # remove extras that may have been queued if we are over the limit
        next
      end


      # Queue nodes after this one
      current_distance = min_distanced_node.distance
      current_positions = min_distanced_node.positions

      # Find neighbouring nodes
      neighbours = nil
      if options[:neighbour_finder]
        neighbours = options[:neighbour_finder].neighbours(min_distanced_node)
      else
        neighbours = min_distanced_node.next_neighbours(graph)
      end

      # explore each neighbour node
      neighbours.each do |onode|
        found_nodes << onode.node.node_id
        length_to_new = 0
        if options[:neighbour_finder]
          # Don't use negative distances as this algorithm cannot handle it, and it is impossible
          # anyway
          if onode.distance > 0
            length_to_new += onode.distance
          else
            length_to_new += 0
          end
        end
        unless first_node
          length_to_new += min_distanced_node.node.length_alone
        end
        new_distance = current_distance + length_to_new
        new_positions = {}
        current_positions.each_pair do |distance, multiplicity|
          new_positions[distance + length_to_new] = multiplicity
        end

        log.debug "Queuing new distance for neighbour: #{onode}: #{new_distance}" if log.debug?
        # new shortest distance found. queue it up
        distanced_node = DistancedOrientedNode.new
        distanced_node.node = onode.node
        distanced_node.first_side = onode.first_side
        distanced_node.distance = new_distance
        distanced_node.positions = new_positions
        counter.push distanced_node

        if options[:ignore_directions]
          log.debug "Queuing new distance for reversed neighbour: #{onode.reverse}: #{new_distance}" if log.debug?
          reverse = DistancedOrientedNode.new
          reverse.node = onode.node
          reverse.first_side = onode.reverse.first_side
          reverse.distance = new_distance
          reverse.positions = new_positions
          counter.push reverse
        end
      end

      first_node = false
    end

    to_return = counter.visited_multiplicity_at_position

    if options[:max_nodes] and found_nodes.length > options[:max_nodes]
      distances_direction_agnostic = {}
      to_return.each do |key, positions|
        prev = distances_direction_agnostic[key[0]]
        distance = positions.keys.min
        if prev.nil? or prev > distance
          distances_direction_agnostic[key[0]] = distance
        end
      end
      if distances_direction_agnostic.length > options[:max_nodes]
        sorted = distances_direction_agnostic.to_a.sort{|a,b| a[1]<=>b[1]}
        # deal with ties i.e. at the edge there can be multiple neighbours
        last_distance = sorted[options[:max_nodes]][1]

        # only keep those nodes that are sufficiently close
        to_return.select! do |key, positions|
          positions.keys.min <= last_distance
        end
      end
    end

    # if ignore directions, then fixup the return so that each direction is included
    if options[:ignore_directions]
      new_to_return = {}
      log.debug "Computing directionless distances for #{to_return}." if log.debug?
      to_return.keys.collect{|key| key[0]}.uniq.each do |node_id|
        keys = [
          Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST,
          Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST].collect do |direction|
            [node_id, direction]
          end
        log.debug "Merging distances for keys #{keys}." if log.debug?
        new_positions = keys.collect{|k| to_return[k]}.reject{|pn| pn.nil?}.reduce do |memo, positions|
          log.debug "Adding directionless distances from #{positions} to #{memo}." if log.debug?
          positions.each_pair do |distance, multiplicity|
            if memo.has_key? distance
              memo[distance] += multiplicity #add ways to get to node at position
            else
              memo[distance] = multiplicity
            end
          end
          memo
        end
        log.debug "Distances merged: #{new_positions}." if log.debug?
        keys.each do |key|
          new_to_return[key] = new_positions
        end
      end
      to_return = new_to_return
    end

    return counter.by_position(to_return)
  end

  # Queue for DistancedOrientedNodes that prioritises nodes by shortest seen distance to node.
  # DistancedOrientedNodes can contain 'position' attribute with a Hash [distance]=>multiplicity which
  # describes a mapping of the start of the node to a position at a set "distance" from an origin node, as
  # well as a "multiplicity" of the mapping: the number of possible paths containing that node at the
  # position.
  # DistancedNodeCounter tracks positions for both visited node, as well as currently queued node,
  # and merges positions of newly pushed nodes. If duplicate nodes are pushed for a current node, the queued
  # node distance and positions are updated when dequeuing.
  # Algorithm is similar to the 'train' algorithm used for the
  # Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble class, with stations counting the number of carriages
  # seen.
  class DistancedNodeCounter
    #include Bio::FinishM::Logging
    attr_reader :visited_multiplicity_at_position

    def initialize
      @visited_multiplicity_at_position = {}
      @current_multiplicity_at_position = {}
      @pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
    end

    def push(dnode)
      # merge position hash for visited node
      visited_positions = @visited_multiplicity_at_position[dnode.to_settable] || {}
      dnode.positions.each_pair do |distance, multiplicity|
        if visited_positions.has_key? distance
          visited_positions[distance] += multiplicity #found m new paths to this position
          #log.debug "Visited node #{dnode.node.node_id} already found at distance #{distance}. Updating multiplicity to #{visited_positions[distance]}." if log.debug?
        else
          visited_positions[distance] = multiplicity
          #log.debug "Visited node #{dnode.node.node_id} found at new distance #{distance}, with multiplicity #{multiplicity}." if log.debug?
        end
      end
      @visited_multiplicity_at_position[dnode.to_settable] = visited_positions

      if @current_multiplicity_at_position.has_key? dnode.to_settable
        #log.debug "Visited node #{dnode.node.node_id} is already queued." if log.debug?
        # node already queued, merge position hash for current node
        current_positions = @current_multiplicity_at_position[dnode.to_settable]
        dnode.positions.each_pair do |distance, multiplicity|
          if current_positions.has_key? distance
            current_positions[distance] += multiplicity
            #log.debug "Queued node already found at distance #{distance}. Updating multiplicity to #{current_positions[distance]}." if log.debug?
          else
            current_positions[distance] = multiplicity
            #log.debug "Queued node found at new distance #{distance}, with multiplicity #{multiplicity}." if log.debug?
          end
        end
      else
        # add position hash and queue node
        @current_multiplicity_at_position[dnode.to_settable] = dnode.positions
        @pqueue.push(dnode, shortest_distance(dnode))
      end
    end

    def shortest_distance(dnode)
      raise unless @visited_multiplicity_at_position.has_key? dnode.to_settable
      @visited_multiplicity_at_position[dnode.to_settable].keys.min
    end

    def shift
      dnode = @pqueue.shift
      unless dnode.nil?
        dnode.positions = @current_multiplicity_at_position.delete dnode.to_settable
        dnode.distance = shortest_distance(dnode) #minimum distance to node
      end
      return dnode
    end

    def by_position(visited_multiplicity_at_position)
      nodes_at_position = {}
      visited_multiplicity_at_position.each_pair do |key, positions|
        positions.each_pair do |distance, multiplicity|
          node_hash = nodes_at_position[distance] || {}
          node_hash[key] = multiplicity
          nodes_at_position[distance] = node_hash
        end
      end
      return nodes_at_position
    end
  end

  # An oriented node some distance from the origin of exploration
  class DistancedOrientedNode
    attr_accessor :node, :first_side, :distance, :positions

    # Using Set object, often we want two separate objects to be considered equal even if
    # they are distinct objects
    def to_settable
      [@node.node_id, @first_side]
    end

    # Which side of the node is not first?
    def second_side
      @first_side == Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST ?
      Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST :
        Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    end

    def next_neighbours(graph)
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
      onode.node = @node
      onode.first_side = @first_side
      return onode.next_neighbours(graph)
    end

    def inspect
      "DistancedOrientedNode #{object_id}: node=#{@node.node_id} first=#{@first_side} distance=#{@distance}"
    end
    alias_method :to_s, :inspect
  end

  class DistancedOrientedTrail
    attr_accessor :oriented_trail, :distance

    def last
      @oriented_trail.last
    end

    def length
      @oriented_trail.length
    end

    def to_settable
      last.to_settable
    end

    def inspect
      @oriented_trail.inspect
    end

    def neighbours_of_last_node(graph)
      last.next_neighbours(graph)
    end
  end
end

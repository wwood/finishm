require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::Dijkstra
  include Bio::FinishM::Logging

  def closest_first_search(graph, oriented_node, options={})
    discovered_oriented_nodes = Set.new

    first = DistancedOrientedNode.new
    first.node = oriented_node.node
    first.first_side = oriented_node.first_side
    first.distance = 0

    pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
    pqueue.push first, first.distance

    first_node = true

    while min_distanced_node = pqueue.shift


      log.debug "Working from #{min_distanced_node.inspect}" if log.debug?
      if discovered_oriented_nodes.include? min_distanced_node.to_settable
        # Already seen this node, do nothing with it
        log.debug "Skipping #{min_distanced_node.inspect} since that has already been seen" if log.debug?
        next
      else
        log.debug "That's a new node, #{min_distanced_node.inspect}" if log.debug?
        # Found a new node for the user to play with
        discovered_oriented_nodes << min_distanced_node.to_settable

        yield min_distanced_node

        #return unless continue

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

          log.debug "Queuing new distance for neighbour: #{onode}: #{new_distance}" if log.debug?
          # new shortest distance found. queue it up
          distanced_node = DistancedOrientedNode.new
          distanced_node.node = onode.node
          distanced_node.first_side = onode.first_side
          distanced_node.distance = new_distance
          pqueue.push distanced_node, new_distance

          if options[:ignore_directions]
            reverse = DistancedOrientedNode.new
            reverse.node = onode.node
            reverse.first_side = onode.reverse.first_side
            reverse.distance = new_distance
            pqueue.push reverse, new_distance
          end
        end

        first_node = false
      end
    end
  end

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

    to_return = {}
    found_nodes = Set.new([initial_oriented_node.node.node_id])
    last_distance = nil

    closest_first_search(graph, initial_oriented_node, options) do |min_distanced_node|

      # Add/overwrite the current one
      to_return[min_distanced_node.to_settable] = min_distanced_node.distance

      #
      found_nodes << min_distanced_node.node.node_id

      if options[:leash_length] and min_distanced_node.distance > options[:leash_length]
        # we are passed leash length, and this is the nearest node. So we are finito.
        log.debug "passed the leash length, cutting short our travels" if log.debug?
        break
      end

      if options[:max_nodes] and found_nodes.length > options[:max_nodes]
        log.debug "passed max-nodes threshold and have #{found_nodes.length} nodes" if log.debug?
        last_distance = min_distanced_node.distance
      end

      if not last_distance.nil? and min_distanced_node.distance > last_distance
        log.debug "found a node further than the max-node threshold" if log.debug?
        break
      end
    end

    if not last_distance.nil?
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

  # An oriented node some distance from the origin of exploration
  class DistancedOrientedNode
    attr_accessor :node, :first_side, :distance

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
end

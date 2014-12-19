require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::Dijkstra
  include Bio::FinishM::Logging

  # Return an array of DistancedOrientedNode objects, those reachable from
  # the initial_oriented_node. options[:leash_length] => max distance explored,
  # can be set to nil to search indefinitely. options[:ignore_directions] =>
  # true or false (default). If true, explore direction-independently.
  # i.e. if 1s->3s and 2s->3s, then include 2s in the returned set of min_distances
  # and continue exploring from 2s.
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

    while min_distanced_node = pqueue.shift

      # Add/overwrite the current one
      to_return[min_distanced_node.to_settable] = min_distanced_node.distance

      log.debug "Working from #{min_distanced_node.inspect}" if log.debug?

      if options[:leash_length] and min_distanced_node.distance > options[:leash_length]
        # we are passed leash length, and this is the nearest node. So we are finito.
        log.debug "passed the leash length, cutting short our travels" if log.debug?
        break
      end

      # Queue nodes after this one
      current_distance = min_distanced_node.distance
      min_distanced_node.next_neighbours(graph).each do |neigh|
        onodes = [neigh]
        onodes.each do |onode|
          new_distance = current_distance
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
      end

      first_node = false
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
      @first_side == OrientedNodeTrail::START_IS_FIRST ?
      OrientedNodeTrail::END_IS_FIRST :
        OrientedNodeTrail::START_IS_FIRST
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
  end
end

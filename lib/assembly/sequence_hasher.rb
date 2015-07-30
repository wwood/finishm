require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::SequenceHasher
  include Bio::FinishM::Logging

  #
  def extend_overlap(graph, oriented_onode, overlap, options={})
    trails = []

    current_path = DistancedOrientedNodeTrail.new
    current_path.add_oriented_node oriented_onode
    current_path.distance = 0

    stack = DS::Stack.new
    stack.push current_path

    # While there is more on the stack
    while current_path = stack.pop

      current_distance = current_path.distance

      if current_distance >= overlap
        # Found all the sequence we need
        trails.push current_path
        next
      end

      # Find neighbouring nodes
      neighbours = nil
      if options[:neighbour_finder]
        neighbours = options[:neighbour_finder].neighbours(oriented_onode)
      else
        neighbours = oriented_node.next_neighbours(graph)
      end

      neighbours.each do |onode|
        new_distance = current_distance
        if options[:neighbour_finder]
          if onode.distance
            new_distance += onode.distance
          else
            new_distance += 0
          end
        end
        new_distance += onode.node.length_alone

        new_path = current_path.copy
        new_path.add_oriented_node onode
        new_path.distance = new_distance
        stack.push new_path
      end
    end
  end

  class DistancedOrientedNodeTrail < Bio::Velvet::Graph::OrientedNodeTrail
    attr_accessor :distance

    def copy
      o = DistancedOrientedNodeTrail.new
      o.trail = Array.new(@trail.collect{|onode| onode.copy})
      o.distance = @distance
      return o
    end

    def to_s
      "DistancedOrientedTrail: #{object_id}: #{to_shorthand} distance=#{@distance}"
    end
  end

end

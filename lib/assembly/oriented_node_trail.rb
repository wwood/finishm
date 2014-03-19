
module Bio
  module Velvet
    class Graph
      # Return an Array of OrientedNode objects corresponding to all the nodes that are next
      # in the graph.
      def neighbours_of(node, first_side)
        neighbour_nodes = nil
        if first_side == OrientedNodeTrail::START_IS_FIRST
          neighbour_nodes = neighbours_off_end(node)
        else
          neighbour_nodes = neighbours_into_start(node)
        end

        neighbours_with_orientation = []
        neighbour_nodes.each do |neighbour|
          arcs = get_arcs_by_node node, neighbour

          # This if statement entered if two nodes are connected twice,
          # in both directions. Remove one direction as it shouldn't be here
          if arcs.length > 1
            if first_side == OrientedNodeTrail::START_IS_FIRST
              arcs = arcs.select do |arc|
                (arc.begin_node_id == node.node_id and arc.begin_node_direction) or
                (arc.end_node_id == node.node_id and !arc.end_node_direction)
              end
            else
              arcs = arcs.select do |arc|
                (arc.end_node_id == node.node_id and arc.end_node_direction) or
                (arc.begin_node_id == node.node_id and !arc.begin_node_direction)
              end
            end
          end

          # Sometimes, but rarely, two nodes will be joined more than once, for whatever reason
          arcs.each do |arc|
            oriented = OrientedNodeTrail::OrientedNode.new
            oriented.node = neighbour
            if arc.begin_node_id == arc.end_node_id
              # A node connecting to itself. Happens rarely.
              if first_side == OrientedNodeTrail::START_IS_FIRST
                if arc.begin_node_direction and arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::START_IS_FIRST
                elsif arc.begin_node_direction and !arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::END_IS_FIRST
                elsif !arc.begin_node_direction and arc.end_node_direction
                  raise "I don't think this is supposed to be possible. Programming error?"
                elsif !arc.begin_node_direction and !arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::START_IS_FIRST
                else
                  raise "programming error"
                end
              else
                # coming from the end of the original node
                if arc.begin_node_direction and arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::END_IS_FIRST
                elsif arc.begin_node_direction and !arc.end_node_direction
                  raise "I don't think this is supposed to be possible. Programming error?"
                elsif !arc.begin_node_direction and arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::START_IS_FIRST
                elsif !arc.begin_node_direction and !arc.end_node_direction
                  oriented.first_side = OrientedNodeTrail::END_IS_FIRST
                else
                  raise "programming error"
                end
              end

            elsif arc.begin_node_id == neighbour.node_id
              # connected to a different node, the 1st in the arc's pair
              if arc.begin_node_direction
                oriented.first_side = OrientedNodeTrail::END_IS_FIRST
              else
                oriented.first_side = OrientedNodeTrail::START_IS_FIRST
              end
            elsif arc.end_node_id == neighbour.node_id
              # connected to a different node, the 2nd in the arc's pair
              if arc.end_node_direction
                oriented.first_side = OrientedNodeTrail::START_IS_FIRST
              else
                oriented.first_side = OrientedNodeTrail::END_IS_FIRST
              end
            end

            neighbours_with_orientation.push oriented
          end
        end

        return neighbours_with_orientation
      end

      # Like #neighbours_of, except takes a velvet node id rather than a node object
      def neighbours_of_node_id(node_id, first_side)
        neighbours_of(@nodes[node_id], first_side)
      end



      # An ordered list of nodes, each with an orientation along that trail
      class OrientedNodeTrail
        include Enumerable
        include Bio::Velvet::Logging

        attr_accessor :trail

        START_IS_FIRST = :start_is_first
        END_IS_FIRST = :end_is_first

        # initialize a new path. If an array is given, each element should be a pair:
        # first element of the pair is a node, and the second true/false or
        # START_IS_FIRST/END_IS_FIRST
        def initialize(node_pairs=[])
          @trail = []
          node_pairs.each do |pair|
            node = pair[0]
            dir = pair[1]
            unless node.kind_of?(Bio::Velvet::Graph::Node) and [true, false, START_IS_FIRST, END_IS_FIRST].include?(dir)
              raise "Bad initialisation of OrientedNodeTrail, with #{node_pairs.inspect}, particularly #{pair.inspect}"
            end
            onode = OrientedNode.new
            onode.node = node
            if dir==true
              onode.first_side = START_IS_FIRST
            elsif dir==false
              onode.first_side = END_IS_FIRST
            else
              onode.first_side = dir
            end
            @trail.push onode
          end
        end

        # Add a node to the trail. start_or_end is either
        # OrientedNodeTrail::START_IS_FIRST or OrientedNodeTrail::END_IS_FIRST
        def add_node(node, start_or_end)
          possible_orientations = [START_IS_FIRST, END_IS_FIRST]
          unless possible_orientations.include?(start_or_end)
            raise "Unexpected orientation in node trail. Need one of #{possible_orientations.inspect}, found #{start_or_end}"
          end
          oriented = OrientedNode.new
          oriented.node = node
          oriented.first_side = start_or_end
          @trail.push oriented
        end

        def add_oriented_node(oriented_node)
          @trail.push oriented_node
        end

        def each(&block)
          @trail.each(&block)
        end

        def last
          @trail[@trail.length-1]
        end

        def remove_last_node
          @trail.pop
        end

        def length
          @trail.length
        end

        def [](index)
          @trail[index]
        end

        # Return true if the path contains the oriented
        # node
        def include_oriented_node?(oriented_node)
          @trail.include?(oriented_node)
        end

        # Return a list of OrientedNode objects, one for each neighbour
        # of the last node in this path (in the correct direction)
        def neighbours_of_last_node(graph)
          graph.neighbours_of(last.node, last.first_side)
        end

        # Return the sequence of the entire trail, or an empty string if there is no
        # nodes in the trail. For certain (small) configurations of (short) nodes, there may
        # be insufficient information to uniquely determine the sequence of the trail.
        # In that case an exception is thrown.
        def sequence
          return '' if @trail.empty?
          twin_nodes_sequence = ''
          fwd_nodes_sequence = ''
          @trail.each do |onode|
            if onode.starts_at_start?
              twin_nodes_sequence = onode.node.ends_of_kmers_of_twin_node + twin_nodes_sequence
              fwd_nodes_sequence += onode.node.ends_of_kmers_of_node
            else
              twin_nodes_sequence = onode.node.ends_of_kmers_of_node + twin_nodes_sequence
              fwd_nodes_sequence += onode.node.ends_of_kmers_of_twin_node
            end
          end
          missing_length_from_each_side = @trail[0].node.parent_graph.hash_length-1
          if twin_nodes_sequence.length < missing_length_from_each_side
            raise Bio::Velvet::NotImplementedException, "Not enough information to know the sequence of a node trail"
          else
            seq_length_required = @trail.collect{|n| n.node.length_alone}.reduce(:+) + missing_length_from_each_side - twin_nodes_sequence.length
            log.debug "first part: #{twin_nodes_sequence}"
            log.debug "second: #{fwd_nodes_sequence[-seq_length_required...fwd_nodes_sequence.length] }"
            return revcom(twin_nodes_sequence)[0...(@trail[0].node.parent_graph.hash_length-1)]+fwd_nodes_sequence
            # calculating this way should be the same, but is somehow buggy in velvet?
            #return revcom(twin_nodes_sequence)+fwd_nodes_sequence[-seq_length_required...fwd_nodes_sequence.length]
          end
        end

        def copy
          o = OrientedNodeTrail.new
          o.trail = Array.new(@trail)
          return o
        end

        def to_s
          "OrientedNodeTrail: #{object_id}: #{collect{|n| [n.node.node_id,n.first_side].join(',')}.join(' ') }"
        end

        def to_short_s
          collect do |onode|
            onode.node.node_id
          end.join(',').to_s
        end

        def inspect
          to_s
        end

        def length_in_bp
          reduce(0){|total, onode| total+=onode.node.length_alone}
        end

        def to_shorthand
          shorthand = @trail.collect do |onode|
            [
              onode.node.node_id,
              onode.starts_at_start? ? 's' : 'e'
              ].join
          end.join(',')
        end

        class OrientedNode
          attr_accessor :node, :first_side

          def initialize(node=nil, first_side=nil)
            @node = node
            if first_side == true
              @first_side = OrientedNodeTrail::START_IS_FIRST
            elsif first_side == false
              @first_side = OrientedNodeTrail::END_IS_FIRST
            else
              @first_side = first_side
            end
          end

          def starts_at_start?
            @first_side == OrientedNodeTrail::START_IS_FIRST
          end

          def starts_at_end?
            @first_side == OrientedNodeTrail::END_IS_FIRST
          end

          def to_s
            "OrientedNode: node #{@node.node_id}, first_side: #{@first_side}"
          end

          # Set#include? doesn't pick up when the same OrientedNode is picked
          # up twice independently, I don't think. So convert to an array first
          def to_settable
            [@node.node_id, @first_side]
          end

          def node_id
            @node.node_id
          end

          def ==(another)
            @node == another.node and @first_side == another.first_side
          end

          def next_neighbours(graph)
            graph.neighbours_of @node, @first_side
          end
        end

        private
        def revcom(seq)
          Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
        end
      end
    end
  end
end

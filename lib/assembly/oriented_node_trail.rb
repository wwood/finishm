module Bio
  module Velvet
    class Graph
      # An ordered list of nodes, each with an orientation along that trail
      class OrientedNodeTrail
        include Enumerable

        START_IS_FIRST = :start_is_first
        END_IS_FIRST = :end_is_first

        def initialize
          @trail = []
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

        def each(&block)
          @trail.each(&block)
        end

        def last
          @trail[@trail.length-1]
        end

        # Return a list of OrientedNode objects, one for each neighbour
        # of the last node in this path (in the correct direction)
        def neighbours_of_last_node(graph)
          neighbour_nodes = nil
          if last.first_side == START_IS_FIRST
            neighbour_nodes = graph.neighbours_off_end(last.node)
          else
            neighbour_nodes = graph.neighbours_into_start(last.node)
          end

          neighbours_with_orientation = neighbour_nodes.collect do |neighbour|
            arcs = graph.get_arcs_by_node last.node, neighbour
            raise "Cycle detected but the code is dumb here currently" if arcs.length != 1 or arcs[0].begin_node_id == arcs[0].end_node_id
            arc = arcs[0]

            oriented = OrientedNode.new
            oriented.node = neighbour
            if arc.begin_node_id == neighbour.node_id
              if arc.begin_node_direction
                oriented.first_side = END_IS_FIRST
              else
                oriented.first_side = START_IS_FIRST
              end
            elsif arc.end_node_id == neighbour.node_id
              if arc.end_node_direction
                oriented.first_side = START_IS_FIRST
              else
                oriented.first_side = END_IS_FIRST
              end
            end

            oriented
          end
          return neighbours_with_orientation
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
            return revcom(twin_nodes_sequence)+fwd_nodes_sequence[-seq_length_required...fwd_nodes_sequence.length]
          end
        end

        def copy
          o = OrientedNodeTrail.new
          each do |oriented|
            o.add_node oriented.node, oriented.first_side
          end
          return o
        end

        def to_s
          "OrientedNodeTrail: #{object_id}: #{collect{|n| [n.node.node_id,n.first_side].join(',')}.join(' ')}"
        end

        def length_in_bp
          reduce(0){|total, onode| total+=onode.node.length_alone}
        end

        class OrientedNode
          attr_accessor :node, :first_side

          def starts_at_start?
            @first_side == OrientedNodeTrail::START_IS_FIRST
          end

          def starts_at_end?
            @first_side == OrientedNodeTrail::END_IS_FIRST
          end

          def to_s
            "OrientedNode: node #{@node.node_id}, first_side: #{@first_side}"
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

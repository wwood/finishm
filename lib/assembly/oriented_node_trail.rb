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

        class OrientedNode
          attr_accessor :node, :first_side

          def starts_at_start?
            @first_side == OrientedNodeTrail::START_IS_FIRST
          end

          def starts_at_end?
            @first_side == OrientedNodeTrail::END_IS_FIRST
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

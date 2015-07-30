
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

        class IllDefinedTrailDefinition < Exception; end
        class InsufficientLengthException < Exception; end

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

        def self.create_from_shorthand(path_string, graph)
          stones = path_string.split(',').collect{|s| s.strip}
          return self.new if stones.length == 0
          trail = []
          stones.each do |stone|
            onode = OrientedNode.new
            if matches = stone.match(/^(\d+)([se])$/)
              node = graph.nodes[matches[1].to_i]
              raise IllDefinedTrailDefinition, "Unable to find node #{matches[1] } in the graph, cannot continue" if node.nil?
              onode.node = node

              if matches[2] == 's'
                onode.first_side = START_IS_FIRST
              else
                onode.first_side = END_IS_FIRST
              end
            else
              raise IllDefinedTrailDefinition, "Unable to underestand shorthand #{stone}"
            end
            trail.push onode
          end
          path = self.new
          path.trail = trail
          return path
        end

        # Given a string like '2,3,4' (super-shorthand form),
        # return the OrientedNodeTrail that thise defines. Raises
        # 'IllDefinedTrailDefinition Exception if there is any ambiguity.
        def self.create_from_super_shorthand(path_string, graph)
          stones = path_string.split(',').collect{|s| s.strip}
          return self.new if stones.length == 0
          if stones.length == 1
            raise IllDefinedTrailDefinition, "Cannot know path orientation when only one node is given"
          end
          state = 'first'
          trail = []

          stones.each do |str|
            if matches = str.match(/^([01-9]+)$/)
              if state == 'first'
                state = 'second'
              elsif state == 'second'
                # Determine the direction of the first two nodes
                first, second = stones[0..1].collect do |str|
                  if matches = str.match(/^([01-9]+)$/)
                    node = graph.nodes[matches[1].to_i]
                    if node.nil?
                      raise IllDefinedTrailDefinition, "Node `#{matches[1] }' from #{path_string} does not appear to be a node ID in the graph"
                    end
                    OrientedNode.new(node, START_IS_FIRST)
                  else
                    raise IllDefinedTrailDefinition, "Unable to parse stepping stone along the path: `#{str}'. Entire path was `#{path_string}'."
                  end
                end
                neighbours_of_first_s = first.next_neighbours(graph)

                rev_first = OrientedNode.new first.node, first.first_side
                rev_first.first_side = END_IS_FIRST
                neighbours_of_first_e = rev_first.next_neighbours(graph)

                if neighbours_of_first_s.find{|n| n.node_id == second.node_id}
                  if neighbours_of_first_e.find{|n| n.node_id == second.node_id}
                    raise IllDefinedTrailDefinition, "Both start and end of first node connect to second node, I'm confused."
                  else
                    seconds = neighbours_of_first_s.select{|n| n.node_id == second.node_id}
                    if seconds.length > 1
                      raise IllDefinedTrailDefinition, "first node connects to both start and end of second node, I'm confused."
                    else
                      trail.push first
                      trail.push seconds[0]
                    end
                  end
                elsif neighbours_of_first_e.find{|n| n.node_id == second.node_id}
                  seconds = neighbours_of_first_e.select{|n| n.node_id == second.node_id}
                  if seconds.length > 1
                    raise IllDefinedTrailDefinition, "first node connects to both start and end of second node, I'm confused."
                  else
                    trail.push rev_first
                    trail.push seconds[0]
                  end
                else
                  raise IllDefinedTrailDefinition, "First and second nodes do not appear to be directly connected"
                end
                state = 'beyond'

              else #we are at the third or later node in the path
                last = trail[-1]
                neighbours_of_last = last.next_neighbours(graph)
                nexts = neighbours_of_last.select{|n| n.node_id == matches[1].to_i}
                if nexts.length == 0
                  raise IllDefinedTrailDefinition, "Nodes #{last} and #{matches[1] } do not appear to be connected"
                elsif nexts.length > 1
                  raise IllDefinedTrailDefinition, "Node #{last} connects to both the start and end of #{matches[1] }, I'm confused"
                else
                  trail.push nexts[0]
                  last = nexts[0]
                end
              end

            else #can't regex the text as shorthand stone or super-shorthand stone
              raise "Unable to parse stepping stone along the path: `#{str}'. Entire path was `#{path_string}'."
            end
          end

          to_return = OrientedNodeTrail.new
          to_return.trail = trail
          return to_return
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

        # Given an Array of [node_id, start_or_end] pairs
        # add these to the trail
        def add_setabled_nodes(setabled_nodes, graph)
          setabled_nodes.each do |pair|
            raise "programming error" if pair.length != 2
            add_node graph.nodes[pair[0]], pair[1]
          end
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

        def delete_at(index)
          @trail.delete_at(index)
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
          fwd_nodes_sequence, twin_nodes_sequence = sequences_within_path
          missing_length_from_each_side = @trail[0].node.parent_graph.hash_length-1
          if twin_nodes_sequence.length < missing_length_from_each_side
            raise InsufficientLengthException, "Not enough information to know the sequence of a node trail"
          else
            seq_length_required = @trail.collect{|n| n.node.length_alone}.reduce(:+) + missing_length_from_each_side - twin_nodes_sequence.length
            log.debug "first part: #{twin_nodes_sequence}"
            log.debug "second: #{fwd_nodes_sequence[-seq_length_required...fwd_nodes_sequence.length] }"
            return revcom(twin_nodes_sequence)[0...(@trail[0].node.parent_graph.hash_length-1)]+fwd_nodes_sequence
            # calculating this way should be the same, but is somehow buggy in velvet?
            #return revcom(twin_nodes_sequence)+fwd_nodes_sequence[-seq_length_required...fwd_nodes_sequence.length]
          end
        end

        def sequences_within_path
          return '', '' if @trail.empty?
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
          return fwd_nodes_sequence, twin_nodes_sequence
        end

        def copy
          o = OrientedNodeTrail.new
          o.trail = Array.new(@trail.collect{|onode| onode.copy})
          return o
        end

        def to_s
          "OrientedNodeTrail: #{object_id}: #{to_shorthand }"
        end

        def to_short_s
          collect do |onode|
            onode.node.node_id
          end.join(',').to_s
        end

        def inspect
          to_s
        end

        # Length of a contig made from this path
        def length_in_bp
          return 0 if @trail.empty?
          return length_in_bp_within_path+@trail[0].node.parent_graph.hash_length-1
        end

        # Length of this trail if it is part of a larger path
        def length_in_bp_within_path
          return 0 if @trail.empty?
          reduce(0) do |total, onode|
            total + onode.node.length_alone
          end
        end

        def to_shorthand
          shorthand = @trail.collect do |onode|
            [
              onode.node.node_id,
              onode.starts_at_start? ? 's' : 'e'
              ].join
          end.join(',')
        end

        def reverse!
          @trail.reverse!
          @trail.each do |onode|
            onode.reverse!
          end
          nil
        end

        def reverse
          rev = copy
          rev.reverse!
          return rev
        end

        # The weighted average of coverages along the trail,
        # (weighted by node length)
        def coverage
          total_length = 0
          total_coverage = 0.0
          each do |onode|
            len =  onode.node.length_alone
            total_coverage += onode.node.coverage*len
            total_length += len
          end
          return total_coverage / total_length
        end

        def ==(another)
          return false if trail.length != another.trail.length
          each_with_index do |onode, i|
            return false unless onode == another[i]
          end
          return true
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

          def to_shorthand
            if @first_side == OrientedNodeTrail::START_IS_FIRST
              return "#{node_id}s"
            else
              return "#{node_id}e"
            end
          end

          # Set#include? doesn't pick up when the same OrientedNode is picked
          # up twice independently, I don't think. So convert to an array first
          def to_settable
            [@node.node_id, @first_side]
          end

          def hash
            to_settable.hash
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

          # switch @first_side of this node
          def reverse!
            if @first_side == OrientedNodeTrail::START_IS_FIRST
              @first_side = OrientedNodeTrail::END_IS_FIRST
            elsif @first_side == OrientedNodeTrail::END_IS_FIRST
              @first_side = OrientedNodeTrail::START_IS_FIRST
            else
              raise "programming error"
            end
          end

          # Return a new OrientedNode with the reverse direction
          def reverse
            rev = OrientedNode.new(@node, @first_side)
            rev.reverse!
            return rev
          end

          def copy
            OrientedNode.new(@node, @first_side)
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

require 'bio-velvet'
require 'bio'

module Bio
  module AssemblyGraphAlgorithms
    class GraphWalkingException < Exception; end

    # A class to extract the sequence of a trail of nodes, in a very lazy, give up easily,
    # kind of way, whenever things get too complicated. This class will throw its hands in
    # the air (by raising and GraphWalkingException) when both ends of the one node connect
    # to an adjacent node on the trail.
    #
    # Hopefully this method should be sufficient for most cases. It might be improved by
    # using a fancier algorithm that considers more than two nodes at a time.
    class LazyGraphWalker
      # Given a list of nodes, one that Node objects when the #each method is called,
      # Return the sequence (as a plain string) of the nodes concatenated together,
      # being cogent of the directionality of the arcs between the nodes.
      #
      # This may not be straightforward, particularly in the presence of palindromic
      # nodes. If any difficulties are encountered, a GraphWalkingException is thrown.
      def trail_sequence(velvet_graph, ordered_collection_of_nodes)
        seq = ''
        last_node = nil
        last_node_was_revcom = nil
        state = :first
        log = Bio::Log::LoggerPlus['finishm']

        ordered_collection_of_nodes.each do |node|
          if last_node.nil?
            # First node in the trail, can't figure it out here
            last_node = node
            state = :second
          else
            # Now there is two nodes in the frame. What is the edge between them, in terms of its direction?
            arcs = velvet_graph.get_arcs_by_node(last_node, node)

            if arcs.empty?
              raise GraphWalkingException, "Attempted to find trail between two unconnected nodes: #{last_node.inspect}, #{node.inspect}"
            elsif arcs.length > 1
              raise GraphWalkingException, "Two adjacent nodes in the graph are (at least) doubly connected too each other, LazyGraphWalker is throwing its hands in the air. Nodes are (in the same order as the specified trail) #{last_node.inspect}, #{node.inspect}"

            # There is a link previous => current, the easy case
            elsif arcs[0].begin_node_id == last_node.node_id
              # There is a connection between the end of node1 and the start of node2,
              # but one or both may be reverse complemented
              arc = arcs[0]
              log.debug "arc last->current found: #{arc.inspect}" if log and log.debug?

              # If the first node has not been written (if this is the second node)
              # then write the sequence of the first node. Here there is a first=>second link,
              # not second=>first, so the forward sequence of the first sequence is written.
              if state == :second
                first_seq = last_node.sequence
                if arc.directions_opposing?
                  first_seq = revcom first_seq
                  log.debug "first node revcom: #{first_seq}" if log and log.debug?
                  seq = first_seq
                else
                  log.debug "first node fwd: #{first_seq}" if log and log.debug?
                  seq = first_seq
                end
                last_node_was_revcom = arc.directions_opposing?
                state = :after_second
              end
              log.debug "At crossroads.\n"+
                "   directions_oppose: #{arc.directions_opposing?}\n"+
                "   last node revcom: #{last_node_was_revcom}" if log and log.debug?

              # Write the sequence of the current node. This will either be the
              # forward sequence or its reverse complement. (Or actually, the abutting bit)
              #
              # If the directions on the arc disagree, then the flow direction is changed
              #
              # The velvet graph structure is setup so that the ends_of_kmers_of_node
              # can be abutted (or ends_of_kmers_of_twin_node) if travelling in the other
              # direction
              #
              # Check that the direction of the current arc and the last arc make sense. Can't have 1=>2, -2=>3, right?
              # Actually you can. That's the same as 1=>2, 2=>-3.
              if arc.directions_opposing?
                if last_node_was_revcom
                  log.
                  seq += node.ends_of_kmers_of_node
                else
                  seq += node.ends_of_kmers_of_twin_node
                end
                last_node_was_revcom = !last_node_was_revcom
              else
                if last_node_was_revcom
                  seq += node.ends_of_kmers_of_twin_node
                else
                  seq += node.ends_of_kmers_of_node
                end
              end
              last_arc = arcs[0]

            # There is a link current => previous
            else
              arc = arcs[0]
              log.debug "arc current->last found: #{arc.inspect}" if log and log.debug?
              if state == :second
                first_seq = last_node.sequence
                if !arc.directions_opposing?
                  first_seq = revcom first_seq
                  log.debug "first node revcom: #{first_seq}" if log and log.debug?
                  seq = first_seq
                else
                  log.debug "first node fwd: #{first_seq}" if log and log.debug?
                  seq = first_seq
                end
                last_node_was_revcom = arc.directions_opposing?
                state = :after_second
              end
              log.debug "At crossroads.\n"+
                "   directions_oppose: #{arc.directions_opposing?}\n"+
                "   last node revcom: #{last_node_was_revcom}" if log and log.debug?

              if arc.directions_opposing?
                if last_node_was_revcom
                  seq += node.ends_of_kmers_of_twin_node
                else
                  seq += node.ends_of_kmers_of_node
                end
                last_node_was_revcom = !last_node_was_revcom
              else
                if last_node_was_revcom
                  seq += node.ends_of_kmers_of_node
                else
                  seq += node.ends_of_kmers_of_twin_node
                end
              end

              last_arc = arcs[0]
            end
            last_node = node
          end
        end
        return '' if last_node.nil? #Return nothing when an empty collection is given.

        if state == :second
          return last_node.sequence
        else
          return seq
        end
      end

      private
      def revcom(seq)
        Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
      end
    end
  end
end
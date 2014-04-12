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
        last_node_used_up_end = nil
        state = :first
        log = Bio::Log::LoggerPlus['finishm']

        add_first_node_forward = lambda do |first_node|
          seq += first_node.sequence
        end
        add_first_node_reverse = lambda do |first_node|
          seq += revcom first_node.sequence
        end

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
            else
              arc = arcs[0]
              # Add the first node, if we are adding the second
              if state == :second
                state = :after_second
                log.debug "Adding start node from arc #{arc}, first_node_id=#{last_node.node_id}, second_node_id=#{node.node_id}" if log and log.debug?
                if arc.connects_end_to_beginning?(last_node.node_id, node.node_id)
                  log.debug "Adding end to beginning" if log and log.debug?
                  add_first_node_forward.call last_node
                  last_node_used_up_end = :start
                elsif arc.connects_end_to_end?(last_node.node_id, node.node_id)
                  log.debug "Adding end to end" if log and log.debug?
                  add_first_node_forward.call last_node
                  last_node_used_up_end = :start
                elsif arc.connects_beginning_to_beginning?(last_node.node_id, node.node_id)
                  log.debug "Adding beginning to beginning" if log and log.debug?
                  add_first_node_reverse.call last_node
                  last_node_used_up_end = :end
                elsif arc.connects_beginning_to_end?(last_node.node_id, node.node_id)
                  log.debug "Adding beginning to end" if log and log.debug?
                  add_first_node_reverse.call last_node
                  last_node_used_up_end = :end
                else
                  raise "Programming error"
                end
              end
              log.debug "At crossroads, with node_id=#{node.node_id}, last_node_used_up_end=#{last_node_used_up_end}, arc=#{arc}" if log and log.debug?

              # Add the new node's sequence
              if state == :after_second
                if arc.connects_end_to_beginning?(last_node.node_id, node.node_id) and last_node_used_up_end == :start
                  log.debug "Adding end to beginning" if log and log.debug?
                  seq += node.ends_of_kmers_of_node
                  last_node_used_up_end = :start
                elsif arc.connects_end_to_end?(last_node.node_id, node.node_id) and last_node_used_up_end == :start
                  log.debug "Adding end to end" if log and log.debug?
                  seq += node.ends_of_kmers_of_twin_node
                  last_node_used_up_end = :end
                elsif arc.connects_beginning_to_beginning?(last_node.node_id, node.node_id) and last_node_used_up_end == :end
                  log.debug "Adding beginning to beginning" if log and log.debug?
                  seq += node.ends_of_kmers_of_node
                  last_node_used_up_end = :start
                elsif arc.connects_beginning_to_end?(last_node.node_id, node.node_id) and last_node_used_up_end == :end
                  log.debug "Adding beginning to end" if log and log.debug?
                  seq += node.ends_of_kmers_of_twin_node
                  last_node_used_up_end = :end
                else
                  raise GraphWalkingException, "The trail being followed to create the trail sequence isn't continuous in a consistent direction. Failed at node #{node.node_id}. last_node_used_up_end=#{last_node_used_up_end}, arc=#{arc}"
                end
              else
                raise "Programming error in state machine"
              end
            end
            last_node = node
          end
        end
        return '' if last_node.nil? #Return nothing when an empty collection is given.

        # When only 1 node is given, return that node's sequence in arbitrary orientation
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

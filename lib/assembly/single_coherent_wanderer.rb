class Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer
  include Bio::FinishM::Logging

  # Like AcyclicConnectionFinder#depth_first_search_with_leash except use
  # single read recoherence. The algorithm used is a generalisation of Dijkstra's
  # shortest path algorithm, where instead of keeping track of the minimum
  # distance to each node, the algorithm keeps track of the distance to a
  # set of nodes long enough to invoke a recoherence kmer.
  def wander(finishm_graph, leash_length, recoherence_kmer, sequence_hash)
    to_return = {}

    # Take the probes and make them all into finishing nodes
    finishing_nodes = []
    finishm_graph.probe_nodes.each_with_index do |probe_node, probe_node_index|
      direction = finishm_graph.probe_node_directions[probe_node_index]
      if direction == true
        finishing_nodes.push [probe_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST]
      else
        finishing_nodes.push [probe_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]
      end
    end

    # Search from each probed node in the graph
    # TODO: is there a better way to implement this by somehow searching with
    # all probe nodes at once, rather than starting fresh with each probe?
    finishm_graph.probe_nodes.each_with_index do |probe_node, probe_node_index|

      # Don't explore from the last node, as no new connections are established
      next if probe_node_index == finishm_graph.probe_nodes.length - 1

      # Go all the way to the leash length,
      # and then search to see if any of the other nodes have been come across
      log.debug "Exploring from probe node \##{probe_node_index+1} (node #{probe_node.node_id}/#{finishm_graph.probe_node_directions[probe_node_index] })" if log.debug?
      pqueue = DS::AnyPriorityQueue.new {|a,b| a < b}
      initial = finishm_graph.initial_path_from_probe(probe_node_index)
      if initial.nil?
        log.warn "Unable to start searching from probe \##{probe_node_index+1}, because it was not found in the graph. Skipping."
        next
      end
      initial_distanced = DistancedOrientedNodeSet.new
      initial_distanced.oriented_trail = initial
      initial_distanced.distance = 0

      # The minimum distance found to get to the head nodes
      minimum_head_nodes_distances = {}
      # Which head node sets is each node connected to?
      node_to_head_node_sets = {}
      #for Logging
      last_logged_node_count = 0

      pqueue.enqueue initial_distanced, 0
      # While there are more node sets in the queue
      while distanced_head_nodes = pqueue.dequeue
        log.debug "Dequeued #{distanced_head_nodes}" if log.debug?
        if log.info? and node_to_head_node_sets.length % 1024 == 0 and node_to_head_node_sets.length > last_logged_node_count
          if last_logged_node_count == 0
            log.info "While exploring from probe \##{probe_node_index+1}.."
          end
          log.info "So far worked with #{node_to_head_node_sets.length} distinct nodes in the assembly graph, at min distance #{distanced_head_nodes.distance}"
          last_logged_node_count = node_to_head_node_sets.length
        end

        settable = distanced_head_nodes.to_settable
        if minimum_head_nodes_distances.key?(settable) and
          distanced_head_nodes.distance >= minimum_head_nodes_distances[distanced_head_nodes.to_settable].distance
          # This node has already been explored, and no shorter path has been found here. Go no further.
          next
        end
        minimum_head_nodes_distances[settable] = distanced_head_nodes
        last_settable = distanced_head_nodes.oriented_trail.last.to_settable
        node_to_head_node_sets[last_settable] ||= Set.new
        node_to_head_node_sets[last_settable] << distanced_head_nodes.to_settable

        if distanced_head_nodes.distance <= leash_length
          # Still within the leash. Push into the stack all the current node's neighbours in the graph
          last = distanced_head_nodes.oriented_trail.last
          neighbour_onodes = finishm_graph.graph.neighbours_of(last.node, last.first_side)
          log.debug "Found #{neighbour_onodes.length} neighbours" if log.debug?
          if neighbour_onodes.length > 1
            # Fork detected. Apply recoherence, and only enqueue those that pass
            log.debug "Multiple neighbours found"
            neighbour_onodes.each do |neighbour|
              candidate = distanced_head_nodes.add_oriented_node_and_copy(neighbour, recoherence_kmer)
              log.debug "Testing recoherence in candidate #{candidate.oriented_trail.to_s}" if log.debug?
              if candidate.last_node_recoherent?(recoherence_kmer, sequence_hash)
                log.debug "Candidate survived recoherence: #{candidate.to_s}" if log.debug?
                pqueue.enqueue candidate, candidate.distance
              elsif log.debug?
                log.debug "Candidate did not survive recoherence #{candidate.oriented_trail.to_s}"
              end
            end
          else
            # One or none neighbours found. Enqueue if there is one
            neighbour_onodes.each do |neighbour|
              candidate = distanced_head_nodes.add_oriented_node_and_copy(neighbour, recoherence_kmer)
              pqueue.enqueue candidate, candidate.distance
            end
          end
        else
          # we are beyond the leash, go no further
        end
      end

      # Now have a hash of minimum distances. Now need to go through those and determine
      # which other nodes the current probe node is connected to
      finishm_graph.probe_nodes.each_with_index do |node, i|
        next if i < probe_node_index # only return the 'upper triangle' of the distance matrices

        finish = finishing_nodes[i]
        heads = node_to_head_node_sets[finish]
        next if heads.nil? #no connection found

        # There might be many head_sets that include the finishing node.
        # Which one has the least distance?
        overall_min_distanced_set = nil
        heads.each do |head_set|
          min_distanced_set = minimum_head_nodes_distances[head_set]
          # If there is a new winner
          if overall_min_distanced_set.nil? or
            overall_min_distanced_set.distance > min_distanced_set.distance

            if probes_on_single_node_ok?(finishm_graph, probe_node_index, i)
              log.debug "Verified that probe indices #{probe_node_index}/#{i} are not failing on a 1 node basis" if log.debug?
            elsif i == probe_node_index+1
              log.debug "Single node probes appear to be circular" if log.debug?
            else
              #TODO: Possibly ok if contigs to be scaffolded are all on the same node. Unlikely in practice due to short tips, but still theoretically possible
              log.debug "Failed to verify that probe indices #{probe_node_index}/#{i} are not failing on a 1 node basis" if log.debug?
              next
            end

            overall_min_distanced_set = min_distanced_set
          end
        end
        next if overall_min_distanced_set.nil? #no connection found - the only connection was a fake one

        min_distance = overall_min_distanced_set.distance
        log.debug "Found a connection between probes #{probe_node_index+1} and #{i+1}, distance: #{min_distance}" if log.debug?
        to_return[[probe_node_index, i]] = min_distance
      end
    end
    return to_return
  end

  # Check for position and orientation if start and finish nodes are
  # on the same velvet node. Return true if OK as below or if the nodes
  # are different
  # -->    <--- OK
  # <--    --> not ok (unless the node is circular)
  # <--    <-- not ok
  # -->    --> not ok
  def probes_on_single_node_ok?(finishm_graph, start_node_index, end_node_index)
    node1 = finishm_graph.probe_nodes[start_node_index]
    node2 = finishm_graph.probe_nodes[end_node_index]
    return true if node1.node_id != node2.node_id

    node1_direction = finishm_graph.probe_node_directions[start_node_index]
    node2_direction = finishm_graph.probe_node_directions[end_node_index]
    node1_offset = direction_independent_offset_of_noded_read_from_start_of_node(
      node1, finishm_graph.probe_node_reads[start_node_index])
    node2_offset = direction_independent_offset_of_noded_read_from_start_of_node(
      node1, finishm_graph.probe_node_reads[end_node_index])
    log.debug "Validating for 1 node problems #{start_node_index}/#{end_node_index} #{node1_direction}/#{node2_direction} offsets #{node1_offset}/#{node2_offset}" if log.debug?

    # true/false and probe1 left of probe2, immediately below, is the most intuitive.
    # but false/true and probe1 right of probe2 is also valid
    if node1_direction == true and node2_direction == false and
      node1_offset < node2_offset
      return true
    end
    if node1_direction == false and node2_direction == true and
      node1_offset > node2_offset
      return true
    end

    if node1_direction == true and node2_direction == false
      onode = finishm_graph.velvet_oriented_node(start_node_index)
      neighbours = finishm_graph.graph.neighbours_of(onode.node, onode.first_side).collect{|n| n.node_id}
      return true if neighbours.include?(node1)
    end

    return false
  end

  private
  def direction_independent_offset_of_noded_read_from_start_of_node(velvet_node, velvet_noded_read)
    if velvet_noded_read.direction == true
      return velvet_noded_read.offset_from_start_of_node
    elsif velvet_noded_read.direction == false
      return velvet_node.corresponding_contig_length - velvet_noded_read.offset_from_start_of_node
    else
      raise "programming error - velvet_noded_read does not have valid direction"
    end
  end

  # An oriented node some distance from the origin of exploration
  class DistancedOrientedNodeSet
    attr_accessor :oriented_trail, :distance

    # Using Set object, often we want two separate objects to be considered equal even if
    # they are distinct objects
    def to_settable
      settable = []
      @oriented_trail.each do |onode|
        settable.push onode.node_id
        settable.push onode.first_side
      end
      return settable
    end

    # Create a copy of this object, then add the given oriented_node
    # to this object, and discard objects from the rear of the trail if they
    # are now of no use for recoherence. Update the distance
    def add_oriented_node_and_copy(oriented_node, recoherence_kmer)
      d = DistancedOrientedNodeSet.new
      new_trail =  @oriented_trail.trail+[oriented_node]

      # Remove unneeded rear nodes that cannot contribute to the recoherence
      # calculation going forward
      cumulative_length = 0
      i = new_trail.length - 1
      while i >= 0 and cumulative_length < recoherence_kmer
        cumulative_length += new_trail[i].node.length_alone
        i -= 1
      end
      i += 1
      d.oriented_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
      d.oriented_trail.trail = new_trail[i..-1]
      # Update distance
      d.distance = @distance+oriented_node.node.length_alone

      return d
    end

    # Is the head nodes single recoherent? Return false if not, otherwise true
    def last_node_recoherent?(recoherence_kmer, sequence_hash)
      @@single_recoherencer ||= Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
      return @@single_recoherencer.validate_last_node_of_path_by_recoherence(
        @oriented_trail,
        recoherence_kmer,
        sequence_hash
        )
    end

    def to_s
      "#{@oriented_trail.to_s}(#{@distance})"
    end
  end
end

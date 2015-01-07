class Bio::FinishM::ProbedGraph
  include Bio::FinishM::Logging
  attr_accessor :probe_nodes, :probe_node_directions, :probe_node_reads, :graph

  attr_accessor :velvet_result_directory

  # Most likely a BinarySequenceStore
  attr_accessor :velvet_sequences

  # Most likely a ReadToNode object
  attr_accessor :read_to_nodes

  # Were all the probe recovered through the process?
  def completely_probed?
    !(@probe_nodes.find{|node| node.nil?})
  end

  def missing_probe_indices
    missings = []
    @probe_nodes.each_with_index do |probe, i|
      missings.push(i+1) if probe.nil?
    end
    return missings
  end

  # Make a Bio::Velvet::Graph::OrientedNodeTrail with just one
  # step in it - the node that corresponds to the probe_index
  def initial_path_from_probe(probe_index)
    initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
    node = @probe_nodes[probe_index]
    raise "No node found for probe #{probe_index}" if node.nil?
    direction = @probe_node_directions[probe_index]

    way = direction ?
    Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
      Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
    initial_path.add_node node, way
    return initial_path
  end

  # Return a Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode
  # corresponding to the index of the probe and its direction
  def velvet_oriented_node(probe_index)
    node = @probe_nodes[probe_index]
    if node.nil?
      return nil
    else
      return initial_path_from_probe(probe_index)[0]
    end
  end

  # The leash is the number of base pairs from the start of the probe,
  # but the path finding algorithm simply uses the combined length of all
  # the nodes without reference to the actual probe sequence. So if the
  # probe is near the end of a long node, then path finding may fail.
  # So adjust the leash length to account for this (or keep the nil
  # if the starting_leash_length is nil)
  def adjusted_leash_length(probe_index, starting_leash_length)
    return nil if starting_leash_length.nil?

    read = @probe_node_reads[probe_index]
    return read.offset_from_start_of_node+starting_leash_length
  end

  # Return a new ProbedGraph that is the same as the current one
  # except that only probe specified in the given probe_indices enumerable
  # are accepted
  def subgraph(probe_indices)
    to_return = Bio::FinishM::ProbedGraph.new
    to_return.graph = @graph
    to_return.velvet_result_directory = @velvet_result_directory
    to_return.velvet_sequences = @velvet_sequences

    to_return.probe_nodes = []
    to_return.probe_node_directions = []
    to_return.probe_node_reads = []
    probe_indices.each do |i|
      to_return.probe_nodes.push @probe_nodes[i-1]
      to_return.probe_node_directions.push @probe_node_directions[i-1]
      to_return.probe_node_reads.push @probe_node_reads[i-1]
    end

    return to_return
  end

  # Return a list of node IDs that are connected through paired-end linkages.
  # This method probably belongs in the Node class except that that is
  # in bio-velvet and yet requires sequence_id_to_node_ids_hash. If all reads are
  # single ended then this method always returns []
  def paired_nodes(node)
    to_return_node_ids = Set.new
    binding.pry
    log.debug "Found #{node.short_reads.length} short reads associated with node #{node}" if log.debug?
    node.short_reads.each do |read|
      pair_read_id = @velvet_sequences.pair_id(read.read_id)
      unless pair_read_id.nil? #i.e. if read is paired
        @read_to_nodes[pair_read_id].each do |node_id|
          to_return_node_ids << node_id
        end
      end
    end
    # Convert node IDs to node objects and return
    return to_return_node_ids.to_a
  end
end

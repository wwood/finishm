class Bio::FinishM::ProbedGraph
  attr_accessor :probe_nodes, :probe_node_directions, :probe_node_reads, :graph

  attr_accessor :velvet_result_directory

  # Most likely a BinarySequenceStore
  attr_accessor :velvet_sequences

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
      to_return.probe_nodes.push @probe_nodes[i]
      to_return.probe_node_directions.push @probe_node_directions[i]
      to_return.probe_node_reads.push @probe_node_reads[i]
    end

    return to_return
  end
end

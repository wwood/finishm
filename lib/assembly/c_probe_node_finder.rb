class Bio::FinishM::CProbeNodeFinder
  include Bio::FinishM::Logging

  # Return an array of [best_node, best_noded_read] that represent the probes in the graph
  def find_probes(velvet_underground_graph, probe_read_ids)
    # First use the C method to extract the set of interesting nodes
    log.debug "Extracting target nodes using the C method.." if log.debug?
    target_node_ids = find_probe_nodes(velvet_underground_graph, probe_read_ids)

    # Then iterate over just those nodes we know are interesting
    log.debug "Extracting from only those #{target_node_ids.length} nodes that are interesting.." if log.debug?
    target_node_ids_set = Set.new target_node_ids
    return Bio::AssemblyGraphAlgorithms::NodeFinder.new.find_unique_nodes_with_sequence_ids(
      velvet_underground_graph, probe_read_ids, :target_node_ids => target_node_ids_set
      )
  end

  # Return a minimal Array of node IDs that contain all probe read IDs
  def find_probe_nodes(velvet_underground_graph, probe_read_ids)
    @bindings ||= Bio::FinishM::VelvetCBinding.new

    c_probe_read_ids = FFI::MemoryPointer.new(:int32, probe_read_ids.length)
    c_probe_read_ids.write_array_of_int32(probe_read_ids)

    probe_nodes = @bindings.extract_best_probe_reads(
      velvet_underground_graph.internal_graph_struct,
      c_probe_read_ids,
      probe_read_ids.length)
    probe_nodes2 = probe_nodes.read_array_of_int32(probe_read_ids.length).collect{|n| n.abs}.uniq

    #clean up
    c_probe_read_ids.free

    return probe_nodes2
  end
end



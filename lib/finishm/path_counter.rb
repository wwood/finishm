class Bio::FinishM::PathCounter
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    options.merge! Bio::FinishM::Visualise::DEFAULT_OPTIONS
    optparse_object.banner = "\nUsage: finishm count_paths --assembly-???

    Count paths through assembly graph
    \n\n"

    optparse_object.separator "Input genome information"
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-exploration arguments:\n\n"
    Bio::FinishM::Visualise.new.add_probe_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options(optparse_object, options)
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    visualise = Bio::FinishM::Visualise.new
    return visualise.validate_argv_length(argv) ||
      visualise.validate_probe_options(options) ||
      visualise.validate_assembly_options(options)
  end

  def run(options, argv)
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    visualise = Bio::FinishM::Visualise.new

    if options[:interesting_probes] or options[:interesting_probe_names]
      finishm_graph, interesting_node_ids, = visualise.generate_graph_from_probes(read_input, options)
      if options[:probe_to_node_map]
        # Output probe map if asked
        visualise.write_probe_to_node_map(options[:probe_to_node_map], finishm_graph, options[:interesting_probes])
      end
    elsif options[:interesting_nodes]
      finishm_graph = visualise.generate_graph_from_nodes(read_input, options)
      interesting_nodes = options[:interesting_nodes]
    elsif options[:assembly_files]
      finishm_graph, interesting_node_ids, = visualise.generate_graph_from_assembly(read_input, options)
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end


    if options[:graph_search_leash_length]
      # get a list of the nodes to be visualised given the leash length
      nodes_within_leash, node_ids_at_leash = visualise.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:graph_search_leash_length] }bp leash"
    else
      nodes_within_leash = finishm_graph.graph.nodes
    end

    initial_onodes = get_leash_start_nodes(finishm_graph, nodes_within_leash)

    log.info "Counting paths through assembly graph.."
    count_paths_through_graph(finishm_graph, initial_onodes,  :range => nodes_within_leash )
  end

  def get_leash_start_nodes(finishm_graph, nodes_within_leash)
    log.info "Finding nodes from which to begin search.."
    start_onodes, = Bio::AssemblyGraphAlgorithms::HeightFinder.new.find_oriented_edge_of_range(finishm_graph.graph, nodes_within_leash)
    if start_onodes.empty?
      #possible with a completely cyclic graph, choose any node to begin
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode finishm_graph.graph.nodes[1], true
      start_onodes = [onode]
      log.info "Graph appears to be a self-contained loop, so chose an arbitrary  node"
    else
      log.info "Found #{start_onodes.length} nodes"
    end
    return start_onodes
  end

  def count_paths_through_graph(finishm_graph, initial_onodes, options={})
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new

    by_height, = height_finder.traverse(finishm_graph.graph, initial_onodes, options)
    min_paths_through = height_finder.min_paths_through(by_height)
    max_paths_through = height_finder.max_paths_through(by_height)
    puts "Minimum number of distinct sequences to explain graph, assuming no errors: #{min_paths_through}."
    puts "Maximum number of distinct sequences allowed by graph: #{max_paths_through}."
  end
end

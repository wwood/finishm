class Bio::FinishM::Visualiser
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm visualise --assembly-??? <output_visualisation_file>

    Visualise an assembly graph
    \n\n"

    options.merge!({
      :graph_search_leash_length => 20000,
      :interesting_probes => nil
    })
    optparse_object.separator "Output visualisation formats (one or more of these must be used)"
    optparse_object.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
      options[:output_graph_png] = arg
    end
    optparse_object.on("--assembly-svg PATH", "Output assembly as a SVG file [default: off]") do |arg|
      options[:output_graph_svg] = arg
    end
    optparse_object.on("--assembly-dot PATH", "Output assembly as a DOT file [default: off]") do |arg|
      options[:output_graph_dot] = arg
    end

    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--probe-ids PROBE_IDS", Array, "explore from these probe IDs in the graph. probe ID is the ID in the velvet Sequence file. See also --leash-length [default: don't start from a node, explore the entire graph]") do |arg|
      options[:interesting_probes] = arg.collect do |read|
        read_id = read.to_i
        if read_id.to_s != read
          raise "Unable to parse probe ID #{read}, from #{arg}, cannot continue"
        end
        read_id
      end
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: unused unless --nodes is specified, otherwise #{options[:graph_search_leash_length]}]") do |arg|
      options[:graph_search_leash_length] = arg
    end

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      if options[:output_graph_png].nil? and options[:output_graph_svg].nil? and options[:output_graph_dot]
        return "No visualisation output format/file given, don't know how to visualise"
      end
      #TODO: this needs to be improved.

      # Need reads unless there is already an assembly
      unless options[:previous_assembly] or options[:previously_serialized_parsed_graph_file]
        return Bio::FinishM::ReadInput.new.validate_options(options, [])
      else
        return nil
      end
    end
  end

  def run(options, argv)
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    # Generate the assembly graph
    log.info "Reading in or generating the assembly graph"
    finishm_graph = nil
    p options[:interesting_probes]
    if options[:interesting_probes]
      dummy_probe_seqs = ['dummy']*options[:interesting_probes].max
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(dummy_probe_seqs, read_input, options)
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end

    # Display the output graph visualisation
    log.info "Converting assembly to a graphviz"
    viser = Bio::Assembly::ABVisualiser.new

    gv = nil
    if options[:interesting_probes].nil?
      gv = viser.graphviz(finishm_graph.graph, {:start_node_ids => finishm_graph.probe_nodes.collect{|node| node.node_id}})
    else
      # Remove nodes unconnected from the interesting read nodes (with appropriate leash length)
      filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
      log.info "Filtering the graph for clarity, removing nodes unconnected to probes #{options[:interesting_probes]}"
      probe_node_ids = options[:interesting_probes].collect do |probe|
        finishm_graph.probe_nodes[probe-1]
      end
      log.debug "ie. filtering velvet nodes unconnected to velvet nodes #{probe_node_ids}"
      filter.remove_unconnected_nodes(finishm_graph.graph,
        probe_node_ids,
        :leash_length => options[:graph_search_leash_length]
        )
      log.info "Filtering finished, leaving #{finishm_graph.graph.nodes.length} nodes and #{finishm_graph.graph.arcs.length} arcs"
      gv = viser.graphviz(finishm_graph.graph, {:start_node_ids => finishm_graph.probe_nodes.collect{|node| node.node_id}})
    end

    if options[:output_graph_png]
      log.info "Writing PNG #{options[:output_graph_png]}"
      gv.output :png => options[:output_graph_png], :use => :neato
    end
    if options[:output_graph_svg]
      log.info "Writing SVG #{options[:output_graph_svg]}"
      gv.output :svg => options[:output_graph_svg], :use => :neato
    end
    if options[:output_graph_dot]
      log.info "Writing DOT #{options[:output_graph_dot]}"
      gv.output :dot => options[:output_graph_dot] if options[:output_graph_dot]
    end
  end
end

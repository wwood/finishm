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
    optparse_object.on("--assembly-svg PATH", "Output assembly as a SVG file [default: off]") do |arg|
      options[:output_graph_svg] = arg
    end
    optparse_object.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
      options[:output_graph_png] = arg
    end
    optparse_object.on("--assembly-dot PATH", "Output assembly as a DOT file [default: off]") do |arg|
      options[:output_graph_dot] = arg
    end

    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--probe-ids PROBE_IDS", Array, "explore from these probe IDs in the graph (comma separated). probe ID is the ID in the velvet Sequence file. See also --leash-length [default: don't start from a node, explore the entire graph]") do |arg|
      options[:interesting_probes] = arg.collect do |read|
        read_id = read.to_i
        if read_id.to_s != read or read_id.nil? or read_id < 1
          raise "Unable to parse probe ID #{read}, from #{arg}, cannot continue"
        end
        read_id
      end
    end
    optparse_object.on("--probe-ids-file PROBE_IDS", String, "explore from the probe IDs given in the file (1 probe ID per line). See also --leash-length [default: don't start from a node, explore the entire graph]") do |arg|
      raise "Cannot specify both --probe-ids and --probe-ids-file sorry" if options[:interesting_probes]
      options[:interesting_probes] = []
      log.info "Reading probe IDs from file: `#{arg}'"
      File.foreach(arg) do |line|
        line.strip!
        next if line == '' or line.nil?
        read_id = line.to_i
        if read_id.to_s != line or read_id < 1 or read_id.nil?
          raise "Unable to parse probe ID #{line}, from file #{arg}, cannot continue"
        end
        options[:interesting_probes].push read_id
      end
      log.info "Read #{options[:interesting_probes].length} probes in"
    end
    optparse_object.on("--node-ids NODE_IDS", Array, "explore from these nodes in the graph (comma separated). Node IDs are the nodes in the velvet graph. See also --leash-length [default: don't start from a node, explore the entire graph]") do |arg|
      options[:interesting_nodes] = arg.collect do |read|
        node_id = read.to_i
        if node_id.to_s != read or node_id.nil? or node_id < 1
          raise "Unable to parse node ID #{read}, from #{arg}, cannot continue"
        end
        node_id
      end
    end
    optparse_object.on("--probe-to-node-map FILE", String, "Output a tab separated file containing the read IDs and their respective node IDs [default: no output]") do |arg|
      options[:probe_to_node_map] = arg
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: unused unless --probe-ids or --nodes is specified, otherwise #{options[:graph_search_leash_length]}]") do |arg|
      options[:graph_search_leash_length] = arg
    end

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      if options[:output_graph_png].nil? and options[:output_graph_svg].nil? and options[:output_graph_dot].nil?
        return "No visualisation output format/file given, don't know how to visualise"
      end
      #TODO: this needs to be improved.
      if options[:interesting_probes] and options[:interesting_nodes]
        return "Can only be interested in probes or nodes, not both, at least currently"
      end

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
    viser = Bio::Assembly::ABVisualiser.new
    gv = nil

    if options[:interesting_probes]
      # Looking based on probes
      if options[:interesting_probes].length > 5
        log.info "Targeting #{options[:interesting_probes].length} probes #{options[:interesting_probes][0..4].join(', ') }, ..."
      else
        log.info "Targeting #{options[:interesting_probes].length} probes #{options[:interesting_probes].inspect}"
      end
      options[:probe_reads] = options[:interesting_probes]
      options[:remove_unconnected_nodes] = true
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

      # Output probe map if asked
      if options[:probe_to_node_map]
        log.info "Writing probe-to-node map to #{options[:probe_to_node_map] }.."
        File.open(options[:probe_to_node_map],'w') do |f|
          f.puts %w(probe_number probe node direction).join("\t")
          finishm_graph.probe_nodes.each_with_index do |node, i|
            if node.nil?
              f.puts [
                i+1,
                options[:interesting_probes][i],
                '-',
                '-',
                ].join("\t")
            else
              f.puts [
                i+1,
                options[:interesting_probes][i],
                node.node_id,
                finishm_graph.probe_node_directions[i] == true ? 'forward' : 'reverse',
                ].join("\t")
            end
          end
        end
      end

      # Create graphviz object
      interesting_node_ids = finishm_graph.probe_nodes.reject{|n| n.nil?}.collect{|node| node.node_id}
      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph, {:start_node_ids => interesting_node_ids})


    elsif options[:interesting_nodes]
      # Looking based on nodes
      if options[:interesting_nodes].length > 5
        log.info "Targeting #{options[:interesting_nodes].length} probes #{options[:interesting_nodes][0..4].join(', ') }, ..."
      else
        log.info "Targeting #{options[:interesting_nodes].length} probes #{options[:interesting_nodes].inspect}"
      end
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

      log.info "Removing unconnected nodes.."
      filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
      interesting_nodes = []
      options[:interesting_nodes].collect do |node_id|
        node  = finishm_graph.graph.nodes[node_id]
        raise "Unable to locate node #{node_id}" if node.nil?
        interesting_nodes.push node
      end
      filter.remove_unconnected_nodes(
        finishm_graph.graph,
        interesting_nodes,
        :leash_length => options[:graph_search_leash_length]
        )

      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph, {:start_node_ids => options[:interesting_nodes]})

    else
      # Visualising the entire graph
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph)
    end

    # Convert gv object to something actually pictorial
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

class Bio::FinishM::Visualise
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm visualise --assembly-??? <output_visualisation_file>

    Visualise an assembly graph
    \n\n"

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

    optparse_object.separator "Input genome information"
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-exploration arguments:\n\n"
    Bio::FinishM::ProbeExplorer.new.add_options(optparse_object, options)

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

      validate_explorer = Bio::FinishM::ProbeExplorer.new.validate_options(options, argv)
      return validate_explorer if validate_explorer

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

    explorer = Bio::FinishM::ProbeExplorer.new
    viser = Bio::Assembly::ABVisualiser.new
    gv = nil

    # Generate the assembly graph
    log.info "Reading in or generating the assembly graph"

    if options[:interesting_probes] or options[:interesting_probe_names]
      # Looking based on probes
      if options[:interesting_probe_names]
        log.info "Targeting #{options[:interesting_probe_names].length} probes through their names e.g. `#{options[:interesting_probe_names] }'"
        options[:probe_read_names] = options[:interesting_probe_names]
      else
        if options[:interesting_probes].length > 5
          log.info "Targeting #{options[:interesting_probes].length} probes #{options[:interesting_probes][0..4].join(', ') }, ..."
        else
          log.info "Targeting #{options[:interesting_probes].length} probes #{options[:interesting_probes].inspect}"
        end
        options[:probe_reads] = options[:interesting_probes]
      end

      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

      # Output probe map if asked
      if options[:probe_to_node_map]
        explorer.write_probe_to_node_map(options[:probe_to_node_map], finishm_graph, options[:interesting_probes])
      end

      # Create graphviz object
      interesting_node_ids = finishm_graph.probe_nodes.reject{|n| n.nil?}.collect{|node| node.node_id}

      nodes_within_leash, node_ids_at_leash = explorer.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:leash_length] }bp leash" if options[:leash_length]

      # Determine paired-end connections
      log.info "Determining paired-end node connections.."
      paired_end_links = explorer.find_paired_end_linkages(finishm_graph, nodes_within_leash)

      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph, {
        :start_node_ids => interesting_node_ids,
        :nodes => nodes_within_leash,
        :end_node_ids => node_ids_at_leash,
        :paired_nodes_hash => paired_end_links,
        })


    elsif options[:interesting_nodes]
      # Looking based on nodes
      if options[:interesting_nodes].length > 5
        log.info "Targeting #{options[:interesting_nodes].length} nodes #{options[:interesting_nodes][0..4].join(', ') }, ..."
      else
        log.info "Targeting #{options[:interesting_nodes].length} node(s) #{options[:interesting_nodes].inspect}"
      end
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

      log.info "Finding nodes within the leash length of #{options[:graph_search_leash_length] }.."
      nodes_within_leash, node_ids_at_leash = explorer.get_nodes_within_leash(finishm_graph, options[:interesting_nodes], options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:leash_length] }bp leash" if options[:leash_length]

      # Determine paired-end connections
      log.info "Determining paired-end node connections.."
      paired_end_links = explorer.find_paired_end_linkages(finishm_graph, nodes_within_leash)

      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph, {
        :start_node_ids => options[:interesting_nodes],
        :nodes => nodes_within_leash,
        :end_node_ids => node_ids_at_leash,
        :paired_nodes_hash => paired_end_links,
        })

    elsif options[:assembly_files]
      # Parse the genome fasta file in
      genomes = Bio::FinishM::InputGenome.parse_genome_fasta_files(
        options[:assembly_files],
        options[:contig_end_length],
        options
        )

      # Create hash of contig end name to probe index
      contig_name_to_probe = {}
      genomes.each do |genome|
        genome.scaffolds.each_with_index do |swaff, scaffold_index|
          probes = [
            genome.first_probe(scaffold_index),
            genome.last_probe(scaffold_index)
            ]
          probes.each do |probe|
            key = nil
            if probe.side == :start
              key = "#{probe.contig.scaffold.name}s"
            elsif probe.side == :end
              key = "#{probe.contig.scaffold.name}e"
            else
              raise "Programming error"
            end

            if contig_name_to_probe.key?(key)
              log.error "Encountered multiple contigs with the same name, this might cause problems, so quitting #{key}"
            end
            contig_name_to_probe[key] = probe.index
          end
        end
      end

      # Gather a list of probe indexes that are of interest to the user
      interesting_probe_ids = []
      if options[:scaffold_sides]
        # If looking at specified ends
        nodes_to_start_from = options[:scaffold_sides].collect do |side|
          if probe = contig_name_to_probe[side]
            interesting_probe_ids << probe
          else
            raise "Unable to find scaffold side in given genome: #{side}"
          end
        end
        log.info "Found #{interesting_probe_ids.length} scaffold sides in the assembly of interest"
      else
        # else looking at all the contig ends in all the genomes
        interesting_probe_ids = contig_name_to_probe.values
        log.info "Visualising all #{interesting_probe_ids.length} contig ends in all genomes"
      end

      # Generate the graph
      probe_sequences = genomes.collect{|genome| genome.probe_sequences}.flatten
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

      # Convert probe IDs into node IDs
      interesting_node_ids = interesting_probe_ids.collect do |pid|
        finishm_graph.probe_nodes[pid].node_id
      end.uniq

      # get a list of the nodes to be visualised given the leash length
      nodes_within_leash, node_ids_at_leash = explorer.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:leash_length] }bp leash" if options[:leash_length]

      # create a nickname hash, id of node to name. Include all nodes even if they weren't specified directly (they only get visualised if they are within leash length of another)
      node_id_to_nickname = {}
      contig_name_to_probe.each do |name, probe|
        key = finishm_graph.probe_nodes[probe].node_id
        if node_id_to_nickname.key?(key)
          node_id_to_nickname[key] += " "+name
        else
          node_id_to_nickname[key] = name
        end
      end

      # Determine paired-end connections
      log.info "Determining paired-end node connections.."
      paired_end_links = explorer.find_paired_end_linkages(finishm_graph, nodes_within_leash)

      # create gv object
      log.info "Converting assembly to a graphviz"
      gv = viser.graphviz(finishm_graph.graph, {
        :start_node_ids => interesting_node_ids,
        :nodes => nodes_within_leash,
        :end_node_ids => node_ids_at_leash,
        :node_id_to_nickname => node_id_to_nickname,
        :paired_nodes_hash => paired_end_links,
        })
    else
      # Visualising the entire graph
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

      # Determine paired-end connections
      log.info "Determining paired-end node connections.."
      # TIM - This will return {} unless get_nodes_within_leash is run first
      paired_end_links = explorer.find_paired_end_linkages(finishm_graph, finishm_graph.graph.nodes)

      log.info "Converting assembly to a graphviz.."
      gv = viser.graphviz(finishm_graph.graph,
        :nodes => finishm_graph.graph.nodes,
        :paired_nodes_hash => paired_end_links)
    end

    # Convert gv object to something actually pictorial
    if options[:output_graph_png]
      log.info "Writing PNG #{options[:output_graph_png] }"
      gv.output :png => options[:output_graph_png], :use => :neato
    end
    if options[:output_graph_svg]
      log.info "Writing SVG #{options[:output_graph_svg] }"
      gv.output :svg => options[:output_graph_svg], :use => :neato
    end
    if options[:output_graph_dot]
      log.info "Writing DOT #{options[:output_graph_dot] }"
      gv.output :dot => options[:output_graph_dot], :use => :neato
    end
  end
end

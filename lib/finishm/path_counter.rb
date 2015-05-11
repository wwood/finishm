class Bio::FinishM::PathCounter
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm count_paths --assembly-???

    Count paths through assembly graph
    \n\n"

    optparse_object.separator "Input genome information"
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-exploration arguments:\n\n"
    Bio::FinishM::ProbeExplorer.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options(optparse_object, options)
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else

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

    # Generate the assembly graph
    log.info "Reading in or generating the assembly graph"

    finishm_graph = nil
    interesting_node_ids = nil
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
      interesting_node_ids = finishm_graph.probe_nodes.reject{|n| n.nil?}.collect{|node| node.node_id}

    elsif options[:interesting_nodes]
      # Looking based on nodes
      if options[:interesting_nodes].length > 5
        log.info "Targeting #{options[:interesting_nodes].length} nodes #{options[:interesting_nodes][0..4].join(', ') }, ..."
      else
        log.info "Targeting #{options[:interesting_nodes].length} node(s) #{options[:interesting_nodes].inspect}"
      end
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
      interesting_node_ids = options[:intereting_nodes]

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

    else
      # Visualising the entire graph
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end

    nodes_within_leash = nil
    if interesting_node_ids
      # get a list of the nodes to be visualised given the leash length
      log.info "Finding nodes within the leash length of #{options[:graph_search_leash_length] }.."
      nodes_within_leash, node_ids_at_leash = explorer.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:leash_length] }bp leash" if options[:leash_length]
    end

    count_paths_through_graph(finishm_graph,{ :range => nodes_within_leash })
  end

  def count_paths_through_graph(finishm_graph, options={})
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, = height_finder.traverse finishm_graph.graph, options
    min_paths_through = height_finder.min_paths_through(by_height)
    max_paths_through = height_finder.max_paths_through(by_height)
    puts "Minimum number of distinct sequences to explain graph, assuming no errors: #{min_paths_through}."
    puts "Maximum number of distinct sequences allowed by graph: #{max_paths_through}."
  end
end

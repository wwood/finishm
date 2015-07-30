class Bio::FinishM::ORFsFinder
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :min_orf_length => 100
    }

  def add_options(optparse_object, options)
    options.merge! Bio::FinishM::Visualise::DEFAULT_OPTIONS
    options.merge! DEFAULT_OPTIONS
    optparse_object.banner = "\nUsage: finishm find_orfs [--orfa OUTPUT_FAA --orfn OUTPUT_FNA]

    Find possible open reading frames in assembly graph
    \n\n"

    optparse_object.separator "Output sequence files"
    add_orfs_options(optparse_object, options)

    optparse_object.separator "Input genome information"
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-exploration arguments:\n\n"
    Bio::FinishM::Visualise.new.add_probe_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options(optparse_object, options)
  end

  def validate_options(options, argv)
    visualise = Bio::FinishM::Visualise.new
    return visualise.validate_argv_length(argv) ||
      visualise.validate_probe_options(options) ||
      visualise.validate_assembly_options(options)
  end

  def add_orfs_options(optparse_object, options)
    optparse_object.on("--orfa OUTPUT_FAA", "Output ORF amino acid sequences [default: off]") do |arg|
      options[:output_faa] = arg
    end
    optparse_object.on("--orfn OUTPUT_FNA", "Output ORF nucleotide sequences [default: orfs.fna unless --orfa is specified]") do |arg|
      options[:output_fna] = arg
    end
  end

  def run(options, argv)
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    visualise = Bio::FinishM::Visualise.new

    if options[:interesting_probes] or options[:interesting_probe_names]
      finishm_graph, interesting_node_ids = visualise.generate_graph_from_probes(read_input, options)
    elsif options[:interesting_nodes]
      finishm_graph = visualise.generate_graph_from_nodes(read_input, options)
      interesting_node_ids = options[:interesting_nodes]
    elsif options[:assembly_files]
      finishm_graph, interesting_node_ids, = visualise.generate_graph_from_assembly(read_input, options)
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end

    if options[:graph_search_leash_length]
      #log.info "Finding nodes within the leash length of #{options[:graph_search_leash_length] }.."
      nodes_within_leash, node_ids_at_leash = visualise.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:graph_search_leash_length] }bp leash" if options[:graph_search_leash_length]

      options[:range] = nodes_within_leash
    else
      options[:range] = finishm_graph.graph.nodes
    end

    initial_onodes = Bio::FinishM::PathCounter.new.get_leash_start_nodes(finishm_graph, options[:range])
    orfs = find_orfs_in_graph(finishm_graph, initial_onodes, options)
    log.info "Found #{orfs.length} open reading frames longer than #{options[:min_orf_length]}."
    if not options[:output_fna] and not options[:output_faa]
        options[:output_fna] = 'orfs.fna'
    end

    if options[:output_fna]
        write_orfs_to_file(orfs, options[:output_fna])
    end
    if options[:output_faa]
        write_orfs_to_file(orfs, options[:output_faa], translate=true)
    end
  end

  def find_orfs_in_graph(finishm_graph, initial_onodes, options={})
    initial_paths = initial_onodes.collect do |onode|
      path = Bio::Velvet::Graph::OrientedNodeTrail.new
      path.add_oriented_node onode
      path
    end

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    orf_trails = orfer.find_orfs_in_graph(finishm_graph.graph, initial_paths,
        options[:min_orf_length], options[:range])

    orfer.orf_sequences_from_trails(orf_trails, options[:min_orf_length])
  end

  def write_orfs_to_file(found_orfs, orfs_file, translate=false)
    if translate
        translator = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    end
    File.open(orfs_file,'w') do |f|
      counter = 0
      found_orfs.each do |name_and_sequence|
        counter += 1
        f.puts ">finishm_orf_#{counter} #{name_and_sequence[0]}"
        if translate
          f.puts translator.sequence2AA(name_and_sequence[1][0...-3])
        else
          f.puts name_and_sequence[1]
        end
      end
    end
  end


  def orf_to_settable(path, start_index, start_offset, end_index, end_offset)
    [path[start_index..end_index].collect{|onode| onode.to_settable},[start_offset, end_offset]]
  end
end

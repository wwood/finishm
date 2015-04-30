class Bio::FinishM::PathCounter
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm pathcounter --assembly-??? --path PATH

    Count paths through assembly graph
    \n\n"

    options.merge!({
    })

    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
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

    count_paths_through_graph(finishm_graph)
  end

  def count_paths_through_graph(finishm_graph, options)
    height_finder = Bio::AssemblyGraphAlgorithms::Height_Finder.new
    by_height, = height_finder.traverse finishm_graph.graph
    min_paths_through = height_finder.min_paths_through(by_height)
    max_paths_through = height_finder.max_paths_through(by_height)
    puts "Minimum number of distinct sequences to explain graph, assuming no errors: #{min_paths_through}."
    puts "Maximum number of distinct sequences allowed by graph: #{max_paths_through}."
  end
end

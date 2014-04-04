class Bio::FinishM::Assembler
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm assemble --assemble-from <node_shorthand> --output-contigs <output.fa> <assembly_definition>

    Assemble
    \n\n"

    options.merge!({
      :output_pathspec => false
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--assemble-from SHORTHAND", "Specify the node and direction to start assembing from e.g. '3s' to start forward from node 3, '4e' to start reverse from node 4 [required]") do |arg|
      unless arg.match(/^\d+[se]+/)
        raise "Unable to parse node shorthand #{arg}"
      end
      options[:initial_node_shorthand] = arg
    end
    optparse_object.on("--output-contigs PATH", "Output found paths to this file in fasta format [required]") do |arg|
      options[:output_trails_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--output-pathspec", "Give the sequence of nodes used in the path in the output contig file [default: #{options[:output_pathspec] }]") do
      options[:output_pathspec] = true
    end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      [
        :initial_node_shorthand,
        :output_trails_file,
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      # Need reads unless there is already an assembly
      unless options[:previous_assembly] or options[:previously_serialized_parsed_graph_file]
        error = Bio::FinishM::ReadInput.new.validate_options(options, [])
        return error unless error.nil?
      else
        return nil
      end
    end
  end

  def run(options, argv)
    assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new

    # Generate the graph
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    graph = finishm_graph.graph

    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand(options[:initial_node_shorthand], graph)
    log.info "Starting to assemble from #{initial_trail.to_shorthand}.."
    path = assembler.assemble_from(initial_trail, graph)

    File.open(options[:output_trails_file],'w') do |output|
      output.print ">#{options[:initial_node_shorthand] }"
      if options[:output_pathspec]
         output.print " #{path.to_shorthand}"
      end
      output.puts
      output.puts path.sequence
    end
  end
end

class Bio::FinishM::Assembler
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm assemble --assemble-from <node_shorthand> --output-contigs <output.fa> <assembly_definition>

    Assemble
    \n\n"

    options.merge!({
      :output_pathspec => false,
      :progressbar => true,
      :min_contig_size => 500,
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--output-contigs PATH", "Output found paths to this file in fasta format [required]") do |arg|
      options[:output_trails_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--assemble-from SHORTHAND", "Specify the node and direction to start assembing from e.g. '3s' to start forward from node 3, '4e' to start reverse from node 4 [default: assemble the entire graph]") do |arg|
      unless arg.match(/^\d+[se]+/)
        raise "Unable to parse node shorthand #{arg}"
      end
      options[:initial_node_shorthand] = arg
    end
    optparse_object.on("--recoherence-kmer LENGTH", Integer, "When paths diverge, try to rescue by using a bigger kmer of this length [default: none]") do |arg|
      options[:recoherence_kmer] = arg
    end
    optparse_object.on("--output-pathspec", "Give the sequence of nodes used in the path in the output contig file [default: #{options[:output_pathspec] }]") do
      options[:output_pathspec] = true
    end
    optparse_object.on("--no-progressbar", "Don't show a progress bar [default: do show one unless --assemble-from is specified]") do
      options[:progressbar] = false
    end
    optparse_object.on("--min-contig-length LENGTH",Integer,"Don't print contigs shorter than this [default: #{options[:min_contig_size] }]") do |arg|
      options[:min_contig_size] = arg
    end
    optparse_object.on("--min-startin-node-coverage COVERAGE",Float,"Only start exploring from nodes with at least this much coverage [default: start from all nodes]") do |arg|
      options[:min_coverage_of_start_nodes] = arg
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
    options[:parse_all_noded_reads] = true if options[:recoherence_kmer]
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    graph = finishm_graph.graph

    if options[:initial_node_shorthand]
      # Just check the shorthand is ok before the time consuming task of reading in more data
      Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand(options[:initial_node_shorthand], graph)
    end

    # Read sequence data
    sequences_file_path = File.join finishm_graph.velvet_result_directory, 'CnyUnifiedSeq'
    log.info "Reading in the actual sequences of all reads from #{sequences_file_path}"
    sequences = Bio::Velvet::Underground::BinarySequenceStore.new sequences_file_path
    log.info "Read in #{sequences.length} sequences"


    if options[:initial_node_shorthand]
      initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand(options[:initial_node_shorthand], graph)
      log.info "Starting to assemble from #{initial_trail.to_shorthand}.."
      path = assembler.assemble_from(initial_trail, graph, sequences,
        :recoherence_kmer => options[:recoherence_kmer]
      )

      File.open(options[:output_trails_file],'w') do |output|
        output.print ">#{options[:initial_node_shorthand] }"
        if options[:output_pathspec]
          output.print " #{path.to_shorthand}"
        end
        output.puts
        output.puts path.sequence
      end
    else

      log.info "Attempting to assemble the entire graph"
      contig_count = 0
      File.open(options[:output_trails_file],'w') do |output|
        progress_io = options[:progressbar] ? $stdout : nil
        assembler.assemble(graph,
          sequences,
          :recoherence_kmer => options[:recoherence_kmer],
          :progressbar_io => progress_io,
          :min_contig_size => options[:min_contig_size],
          :min_coverage_of_start_nodes => options[:min_coverage_of_start_nodes],
          ) do |path|
            contig_count += 1
            output.print ">contig#{contig_count}"
            if options[:output_pathspec]
              output.print " #{path.to_shorthand}"
            end
            output.puts
            output.puts path.sequence
          end
      end
      log.info "Assembled #{contig_count} contigs"
    end
  end
end

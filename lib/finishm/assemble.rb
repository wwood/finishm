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
      :bubbly => false,
      :max_tip_length => Bio::AssemblyGraphAlgorithms::BubblyAssembler::DEFAULT_MAX_TIP_LENGTH,
      :max_bubble_length => Bio::AssemblyGraphAlgorithms::BubblyAssembler::DEFAULT_MAX_BUBBLE_LENGTH,
      :bubble_node_count_limit => Bio::AssemblyGraphAlgorithms::BubblyAssembler::DEFAULT_BUBBLE_NODE_COUNT_LIMIT,
      :min_confirming_recoherence_kmer_reads => Bio::AssemblyGraphAlgorithms::SingleEndedAssembler::DEFAULT_MIN_CONFIRMING_RECOHERENCE_READS,
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
    optparse_object.on("--recoherence-min-reads NUM", Integer, "Number of reads required to agree with recoherence [default: #{options[:min_confirming_recoherence_kmer_reads] } (when --recoherence-kmer is specified)]") do |arg|
      options[:min_confirming_recoherence_kmer_reads] = arg
    end
    optparse_object.on("--max-tip-length LENGTH", Integer, "Maximum length of 'tip' in assembly graph to ignore [default: #{options[:max_tip_length] }]") do |arg|
      options[:max_tip_length] = arg
    end
    optparse_object.on("--bubbly", "Assemble with the bubbly method [default: #{options[:bubbly] }]") do
      options[:bubbly] = true
    end
    optparse_object.on("--max-bubble-size NUM", Integer, "Max bubble size available for bubbly method [default: #{options[:max_bubble_length] }]") do |arg|
      options[:max_bubble_length] = arg
    end
    optparse_object.on("--max-bubble-complexity NUM", Integer, "Max number of nodes in a bubble to explore before giving up (0 for infinite) [default: #{options[:bubble_node_count_limit] }]") do |arg|
      if arg == 0
        options[:bubble_node_count_limit] = nil
      else
        options[:bubble_node_count_limit] = arg
      end
    end
    optparse_object.on("--output-pathspec", "Give the sequence of nodes used in the path in the output contig file [default: #{options[:output_pathspec] }]") do
      options[:output_pathspec] = true
    end
    optparse_object.on("--output-contig-stats FILE", "Output stats about each contig to this file [default: don't output anything]") do |arg|
      options[:output_stats] = arg
    end
    optparse_object.on("--no-progressbar", "Don't show a progress bar [default: do show one unless --assemble-from is specified]") do
      options[:progressbar] = false
    end
    optparse_object.on("--min-contig-length LENGTH",Integer,"Don't print contigs shorter than this [default: #{options[:min_contig_size] }]") do |arg|
      options[:min_contig_size] = arg
    end
    optparse_object.on("--min-starting-node-coverage COVERAGE",Float,"Only start exploring from nodes with at least this much coverage [default: start from all nodes]") do |arg|
      options[:min_coverage_of_start_nodes] = arg
    end
    optparse_object.on("--min-starting-node-length LENGTH",Integer,"Only start exploring from nodes with at least this length [default: start from all nodes]") do |arg|
      options[:min_length_of_start_nodes] = arg
    end
    optparse_object.on("--max-coverage-at-fork COVERAGE",Float,"When reached a fork, don't take paths with more than this much coverage [default: not applied]") do |arg|
      options[:max_coverage_at_fork] = arg
    end
    optparse_object.on("--badformat FILE", "Output contigs in badformat file") do |arg|
      options[:output_badformat_file] = arg
    end
    optparse_object.on("--debug", "Build the graph, then drop to a pry console. [default: #{options[:debug] }]") do
      options[:debug] = true
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
    # Generate the graph
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    finishm_graph = nil
    if options[:recoherence_kmer].nil?
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options.merge({
        :dont_parse_reads => true,
        :dont_parse_noded_reads => true,
        }))
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end
    graph = finishm_graph.graph

    if options[:initial_node_shorthand]
      # Just check the shorthand is ok before the time consuming task of reading in more data
      Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand(options[:initial_node_shorthand], graph)
    end

    # Setup assembler
    assembler = nil
    if options[:bubbly]
      assembler = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
    else
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
    end
    [
      :recoherence_kmer,
      :min_confirming_recoherence_kmer_reads,
      :min_contig_size,
      :min_coverage_of_start_nodes,
      :min_length_of_start_nodes,
      :max_tip_length,
      :leash_length,
      :max_bubble_length,
      :bubble_node_count_limit,
      :max_coverage_at_fork,
      ].each do |opt|
        assembler.assembly_options[opt] = options[opt]
      end
    assembler.assembly_options[:sequences] = finishm_graph.velvet_sequences

    binding.pry if options[:debug]

    if options[:initial_node_shorthand]
      initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand(options[:initial_node_shorthand], graph)
      log.info "Starting to assemble from #{initial_trail.to_shorthand}.."
      path, visited_nodes = assembler.assemble_from(initial_trail)

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
      stats_output = nil
      if options[:output_stats]
        stats_output = File.open(options[:output_stats],'w')
        stats_output.puts %w(name coverage).join("\t")
      end
      badformat_writer = Bio::FinishM::BadFormatWriter.new
      File.open(options[:output_trails_file],'w') do |output|
        progress_io = options[:progressbar] ? $stdout : nil
        assembler.assembly_options[:progressbar_io] = progress_io
        assembler.assemble do |path|
          contig_count += 1
          name = "contig#{contig_count}"
          output.print ">#{name}"
          if options[:output_pathspec]
            output.print " #{path.to_shorthand}"
          end
          output.puts
          output.puts path.sequence

          if !stats_output.nil?
            stats_output.puts [
              name,
              path.coverage,
              ].join("\t")
          end

          badformat_writer.add_metapath(name, path) if options[:output_badformat_file]
        end
      end
      if options[:output_badformat_file]
        File.open(options[:output_badformat_file],'w') do |out|
          badformat_writer.write out
        end
      end
      log.info "Assembled #{contig_count} contigs"
      stats_output.close if !stats_output.nil?
    end
  end
end

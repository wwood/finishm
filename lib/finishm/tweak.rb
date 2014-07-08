class Bio::FinishM::Tweaker
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :contig_end_length => 200,
    :graph_search_leash_length => 20000,
    :unscaffold_first => false,
    :recoherence_kmer => 1,
    }

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm roundup --genomes <genome1.fasta>[,<genome2.fasta>,...] --fastq-gz <reads..> --output-directory <directory>

    Takes one or more genomes and tries to improve their quality by reducing the number of
    scaffolds and N characters they contain.

    Example:

      finishm roundup --genomes genome1.fasta,genome2.fasta --fastq-gz reads.1.fq.gz,reads.2.fq.gz --output-directory finishm_roundup_results

    That will create a collapsed de-Bruijn graph from reads.1.fq.gz and reads.2.fq.gz, then try to find connections between
    the starts and the ends of the contigs in genome1.fasta. If any connections between contigs are mutually exclusive,
    then they are incorporated into scaffolds together, and gapfilling is attempted. The final sequences are output in
    the finishm_roundup_results directory in FASTA format. The procedure is then repeated for genome2.fasta.

    \n\n"

    options.merge!(DEFAULT_OPTIONS)

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--genomes FASTA_1[,FASTA_2...]", Array, "fasta files of genomes to be improved [required]") do |arg|
      options[:assembly_files] = arg
    end
    optparse_object.on("--output-directory PATH", "Output results to this directory [required]") do |arg|
      options[:output_directory] = arg
    end

    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length] }]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--recoherence-kmer NUM", Integer, "Use a kmer longer than the original velvet one, to help remove bubbles and circular paths [default: none]") do |arg|
      options[:recoherence_kmer] = arg
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: #{options[:graph_search_leash_length] }]") do |arg|
      options[:graph_search_leash_length] = arg
    end
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs. This option is only relevant to the wander step; gapfilling is attempted on all gaps by default. [default: #{options[:unscaffold_first] }]") do
      options[:unscaffold_first] = true
    end
    #optparse_object.on("--proceed-on-short-contigs", "By default, when overly short contigs are encountered, finishm croaks. This option stops the croaking [default: #{options[:proceed_on_short_contigs] }]") do
    #  options[:proceed_on_short_contigs] = true
    #end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      [
        :assembly_files,
        :output_directory,
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv=[])
    # Make sure output directory is writeable to avoid late croaking
    output_directory = setup_output_directory options[:output_directory]

    # Gather the probes from each genome supplied
    genomes = Bio::FinishM::InputGenome.parse_genome_fasta_files(
      options[:assembly_files],
      options[:contig_end_length],
      options
      )

    # Generate one velvet assembly to rule them all (forging the assembly is hard work..)
    probe_sequences = genomes.collect{|genome| genome.probe_sequences}.flatten
    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    master_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    gapfill_scaffold = lambda do |gapfiller, genome, scaffold_index|
      connections = []
      genome.each_gap_probe_pair(scaffold_index) do |probe1, probe2|
        log.debug "Gapfilling between probes #{probe1.number} and #{probe2.number}.."
        connections.push gapfiller.gapfill(master_graph, probe1.index, probe2.index, options)
      end
      connections
    end

    revcom = lambda do |seq|
      Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
    end

    # For each genome, wander, gapfill, then output
    wanderer = Bio::FinishM::Wanderer.new
    gapfiller = Bio::FinishM::GapFiller.new
    printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
    genomes.each do |genome|
      # wander using just the probes on the ends of the scaffolds
      connected_scaffolds = wander_a_genome(wanderer, genome, master_graph, options)
      output_path = File.join(output_directory, File.basename(genome.filename)+".scaffolds.fasta")

      File.open(output_path, 'w') do |output_file|
        # gapfill between
        # (1) interpreted_connections
        # (2) gaps that were present before above wander
        connected_scaffolds.each_with_index do |cross_scaffold_connection, connected_scaffold_index|
          pretend_contig = cross_scaffold_connection.contigs[0]
          first_scaffold_index = pretend_contig.sequence_index
          first_scaffold = genome.scaffolds[first_scaffold_index]
          scaffold_sequence = first_scaffold.contigs[0].sequence

          # Gapfill contigs within the scaffold on the extreme LHS
          connections = gapfill_scaffold.call(gapfiller, genome, first_scaffold_index)
          connections.each_with_index do |aconn, i|
            scaffold_sequence = printer.one_connection_between_two_contigs(
              master_graph.graph, scaffold_sequence, aconn, first_scaffold.contigs[i+1].sequence
              )
          end
          scaffold_sequence = revcom.call(scaffold_sequence) if pretend_contig.direction == false

          last_contig = nil
          cross_scaffold_connection.contigs.each_with_index do |contig, superscaffold_index|
            unless last_contig.nil? #skip the first contig - it be done

              # Ready the contig on the RHS of this join
              # Gapfill within the scaffold on the RHS of the new gap
              rhs_sequence = genome.scaffolds[contig.sequence_index].contigs[0].sequence
              second_scaffold_index = cross_scaffold_connection.contigs[superscaffold_index].sequence_index
              connections = gapfill_scaffold.call(gapfiller, genome, second_scaffold_index)
              connections.each_with_index do |aconn, scaffold_index|
                second_sequence = genome.scaffolds[cross_scaffold_connection.contigs[superscaffold_index].sequence_index].contigs[scaffold_index+1].sequence
                rhs_sequence = printer.one_connection_between_two_contigs(
                  master_graph.graph,
                  rhs_sequence,
                  aconn,
                  second_sequence
                  )
              end

              # Gapfill across the new gap between scaffolds
              aconn = gapfiller.gapfill(master_graph,
                last_contig.direction == true ? genome.last_probe(last_contig.sequence_index).index : genome.first_probe(last_contig.sequence_index).index,
                contig.direction == true ? genome.first_probe(contig.sequence_index).index : genome.last_probe(contig.sequence_index).index,
                options
                )
              second_sequence = genome.scaffolds[contig.sequence_index].contigs[0].sequence
              rhs_sequence = revcom.call(rhs_sequence) if contig.direction == false
              scaffold_sequence = printer.one_connection_between_two_contigs(
                master_graph.graph,
                scaffold_sequence,
                aconn,
                rhs_sequence
                )
            end
            last_contig = contig
          end

          #Output the scaffold to the output directory
          scaffold_names = cross_scaffold_connection.contigs.collect do |contig|
            genome.scaffolds[contig.sequence_index].name
          end
          output_file.puts ">scaffold#{connected_scaffold_index+1} #{scaffold_names.join(':') }"
          output_file.puts scaffold_sequence
        end

        log.info "Wrote #{connected_scaffolds.length} scaffolds to #{output_path}"
      end
    end
  end

  def setup_output_directory(given_directory)
    output_directory = File.absolute_path(given_directory)
    log.debug "Using output directory: #{output_directory}" if log.debug?

    if File.exist?(output_directory)
      if !File.directory?(output_directory)
        log.error "Specified --output-directory #{output_directory} exists but is a file and not a directory. Cannot continue."
        exit 1
      elsif !File.writable?(output_directory)
        log.error "Specified --output-directory #{output_directory} is not writeable. Cannot continue."
        exit 1
      else
        log.debug "Already existing output directory #{output_directory} seems usable"
      end
    else
      # Creating a new output directory
      Dir.mkdir(output_directory)
    end

    return output_directory
  end

  def wander_a_genome(wanderer, genome, master_probed_graph, options)
    # Create new finishm_graph with only probes from the ends of the scaffolds of this genome
    probe_indices = []
    genome.each_scaffold_end_numbered_probe{|probe| probe_indices.push(probe.number)}
    genome_graph = master_probed_graph.subgraph(probe_indices)

    num_scaffolds = genome.scaffolds.length

    all_connections = wanderer.probed_graph_to_connections(genome_graph, options)

    interpreter = Bio::FinishM::ConnectionInterpreter.new(all_connections, (0...num_scaffolds))
    connections = interpreter.doubly_single_contig_connections
    log.debug "Found #{connections.length} connections between contigs that can be used for scaffolding" if log.debug?

    return interpreter.scaffolds(connections)
  end
end

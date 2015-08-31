class Bio::FinishM::RoundUp
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :contig_end_length => 200,
    :graph_search_leash_length => 20000,
    :unscaffold_first => false,
    :recoherence_kmer => 1,
    :debug => false,
    :gapfill_only => false,
    :max_explore_nodes => 10000,
    :max_gapfill_paths => 10,
    :gapfill_with_max_coverage => false,
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
    optparse_object.on("--gapfill-only", "Don't wander, just gapfill [default: #{options[:gapfill_only] }]") do
      options[:gapfill_only] = true
    end
    optparse_object.on("--max-gapfill-paths NUM", Integer, "When this number of paths is exceeded, don't gapfill, print as Ns [default: #{options[:max_gapfill_paths] }]") do |arg|
      options[:max_gapfill_paths] = arg
    end
    optparse_object.on("--max-explore-nodes NUM", Integer, "Only explore this many nodes. If max is reached, do not make connections. [default: #{options[:max_explore_nodes] }]") do |arg|
      options[:max_explore_nodes] = arg
    end
    optparse_object.on("--gapfill-with-max-coverage", "When gapfilling, take the path with maximal coverage and do not print variants [default: #{options[:gapfill_with_max_coverage] }]") do
      options[:gapfill_with_max_coverage] = true
    end
    optparse_object.on("--debug", "Build the graph, then drop to a pry console. [default: #{options[:debug] }]") do
      options[:debug] = true
    end
    optparse_object.on("--probe NUM",Integer,"debug mode: explore from this probe number (1-based index), gapfill only, no wander. [default: explore from all probes}]") do |arg|
      options[:interesting_probes] = [arg-1] #convert to 0-based index
      options[:gapfill_only] = true
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

    binding.pry if options[:debug]

    # For each genome, wander, gapfill, then output
    wanderer = Bio::FinishM::Wanderer.new
    gapfiller = Bio::FinishM::GapFiller.new
    printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
    genomes.each do |genome|
      # wander using just the probes on the ends of the scaffolds
      connected_scaffolds = nil
      all_connections = []
      gaps_filled_in_genome = 0
      wandered_probe_indices = nil

      File.open(File.join(output_directory, File.basename(genome.filename)+".report.txt"),'w') do |report|
        report.puts "#{Time.now} FinishM report for roundup run with: #{options.inspect}"

        if options[:gapfill_only]
          log.info "Skipping wander, gapfilling only"
          connected_scaffolds = Bio::FinishM::ConnectionInterpreter.new([], (0...genome.scaffolds.length)).scaffolds([])
        else
          log.debug "Wandering.."
          connected_scaffolds, all_connections, wandered_probe_indices = wander_a_genome(wanderer, genome, master_graph, options, report)
        end

        # Write out all the connections
        File.open(File.join(output_directory, File.basename(genome.filename)+".connections.csv"),'w') do |con_file|
          all_connections.each do |connection|
            con_file.puts connection
          end
        end

        output_path = File.join(output_directory, File.basename(genome.filename)+".scaffolds.fasta")
        variants_path = File.join(output_directory, File.basename(genome.filename)+".at_least_half_completely_wrong.vcf")
        num_circular_scaffolds = 0

        File.open(output_path, 'w') do |output_file|
          File.open(variants_path,'w') do |variants_file|
            variants_file.puts %w(#CHROM POS ID REF ALT QUAL FILTER INFO).join("\t")
            # gapfill between
            # (1) interpreted_connections
            # (2) gaps that were present before above wander
            connected_scaffolds.each_with_index do |cross_scaffold_connection, connected_scaffold_index|
              superscaffold_name = "scaffold#{connected_scaffold_index+1}"

              pretend_contig = cross_scaffold_connection.contigs[0]
              first_scaffold_index = pretend_contig.sequence_index

              # Gapfill contigs within the scaffold on the extreme LHS
              scaffold_sequence, num_gaps, variants = gapfill_a_scaffold(gapfiller, printer, master_graph, genome, first_scaffold_index, pretend_contig.direction, superscaffold_name, report, variants_file, options)
              gaps_filled_in_genome += num_gaps

              last_contig = nil
              cross_scaffold_connection.contigs.each_with_index do |contig, superscaffold_index|
                unless last_contig.nil? #skip the first contig - it be done
                  last_name = genome.scaffolds[last_contig.sequence_index].name
                  current_name = genome.scaffolds[contig.sequence_index].name
                  log.debug "Connecting #{last_name} and #{current_name}" if log.debug?

                  # Ready the contig on the RHS of this join
                  # Gapfill within the scaffold on the RHS of the new gap
                  rhs_sequence, num_gaps, variants = gapfill_a_scaffold(gapfiller, printer, master_graph, genome, contig.sequence_index, contig.direction, superscaffold_name, report, variants_file, options)
                  gaps_filled_in_genome += num_gaps

                  # Gapfill across the new gap between scaffolds
                  aconn = gapfiller.gapfill(master_graph,
                    last_contig.direction == true ? genome.last_probe(last_contig.sequence_index).index : genome.first_probe(last_contig.sequence_index).index,
                    contig.direction == true ? genome.first_probe(contig.sequence_index).index : genome.last_probe(contig.sequence_index).index,
                    options
                    )
                  second_sequence = genome.scaffolds[contig.sequence_index].contigs[0].sequence
                  log.debug "Found #{aconn.paths.length} connections between #{last_name} and #{current_name}" if log.debug?
                  if aconn.paths.length == 0
                    # when this occurs, it is due to there being a circuit in the path, so no paths are printed.
                    # (at least for now) TODO: this could be improved.
                    # Just arbitrarily put in 100 N characters, to denote a join, but no gapfill
                    scaffold_sequence = scaffold_sequence+('N'*100)+rhs_sequence
                  else
                    acon.collapse_paths_to_maximal_coverage_path! if options[:gapfill_with_max_coverage]
                    scaffold_sequence, variants = printer.ready_two_contigs_and_connections(
                      master_graph.graph,
                      scaffold_sequence,
                      aconn,
                      rhs_sequence,
                      master_graph.velvet_sequences
                      )
                    # Print variants
                    # TODO: need to change coordinates of variants, particularly when >2 contigs are joined?
                    variants.each do |variant|
                      variant.reference_name = superscaffold_name
                      variants_file.puts variant.vcf(scaffold_sequence)
                    end
                  end
                end
                last_contig = contig
              end

              #Output the scaffold to the output directory
              descriptor = nil
              if cross_scaffold_connection.circular?
                descriptor = 'circular'
                num_circular_scaffolds += 1
              else
                descriptor = 'scaffold'
              end
              scaffold_names = cross_scaffold_connection.contigs.collect do |contig|
                genome.scaffolds[contig.sequence_index].name
              end
              output_file.puts ">#{superscaffold_name} #{descriptor} #{scaffold_names.join(':') }"
              output_file.puts scaffold_sequence
            end

            num_connected_scaffolds = genome.scaffolds.length - connected_scaffolds.length
            log.info "Wrote #{connected_scaffolds.length} scaffolds to #{output_path}, after scaffolding #{num_connected_scaffolds} scaffolds together (#{num_circular_scaffolds} were circular) and filling #{gaps_filled_in_genome} gaps."
          end
        end
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

  def wander_a_genome(wanderer, genome, master_probed_graph, options, report)
    # Create new finishm_graph with only probes from the ends of the scaffolds of this genome
    probe_indices = []
    genome.each_scaffold_end_numbered_probe{|probe| probe_indices.push(probe.number)}
    genome_graph = master_probed_graph.subgraph(probe_indices)

    num_scaffolds = genome.scaffolds.length

    all_connections = wanderer.probed_graph_to_connections(genome_graph, options)

    interpreter = Bio::FinishM::ConnectionInterpreter.new(all_connections, (0...num_scaffolds))
    connections = interpreter.doubly_single_contig_connections
    report.puts "Found #{connections.length} connections between contigs that can be used for scaffolding"
    unconnected_probes = interpreter.unconnected_probes
    report.puts "Found #{unconnected_probes.length} contig ends that did not connect to any others"
    unconnected_probes.each do |probe|
      report.puts "Did not connect to any other probes: #{probe.inspect}"
    end

    return interpreter.scaffolds(connections), all_connections, probe_indices
  end

  def gapfill_a_scaffold(gapfiller, printer, master_graph, genome, scaffold_index, scaffold_direction, superscaffold_name, report, variants_file, options)
    connections = []
    genome.each_gap_probe_pair(scaffold_index) do |probe1, probe2|
      log.info "Gapfilling between probes #{probe1.number+1} and #{probe2.number+1}.."
      next unless options[:interesting_probes].nil? or
      options[:interesting_probes].include?(probe1.number) or
      options[:interesting_probes].include?(probe2.number)
      connections.push gapfiller.gapfill(master_graph, probe1.index, probe2.index, options)
    end
    log.debug "Found #{connections.length} connections" if log.debug?

    all_variants = []
    num_gapfills = 0
    scaffold = genome.scaffolds[scaffold_index]
    gapfilled_sequence = genome.scaffolds[scaffold_index].contigs[0].sequence
    connections.each_with_index do |aconn, i|
      rhs_sequence = scaffold.contigs[i+1].sequence
      gapfilled_sequence, variants, gapfilled = piece_together_gapfill(
        printer, master_graph, gapfilled_sequence, aconn, rhs_sequence, genome.gap_length(scaffold_index, i),
        options[:max_gapfill_paths]
        )
      if gapfilled
        num_gapfills += 1
        variants.each{|v| all_variants << v}
        to_log = "Filled a gap on genome #{genome.filename}: scaffold #{scaffold.name}: #{scaffold.contigs[i].scaffold_position_end+1}-#{scaffold.contigs[i+1].scaffold_position_start-1}"
        report.puts to_log
        log.info to_log
      end
    end
    if scaffold_direction == false
      gapfilled_sequence = revcom(gapfilled_sequence)
      all_variants.each do |variant|
        variant.position = gapfilled_sequence.length - variant.position
        variant.reverse!
      end
    end
    all_variants.each do |variant|
      variant.reference_name = superscaffold_name
      variants_file.puts variant.vcf(gapfilled_sequence)
    end
    return gapfilled_sequence, num_gapfills, all_variants
  end

  def piece_together_gapfill(printer, master_graph, first_sequence, aconn, second_sequence, gap_length, max_gapfill_paths, options={})
    scaffold_sequence = nil
    gapfilled = -1
    if aconn.paths.length == 0 or aconn.paths.length > max_gapfill_paths
      # No paths found. Just fill with Ns like it was before
      scaffold_sequence = first_sequence + 'N'*gap_length + second_sequence
      gapfilled = false
    else
      acon.collapse_paths_to_maximal_coverage_path! if options[:gapfill_with_max_coverage]
      scaffold_sequence, variants = printer.ready_two_contigs_and_connections(
        master_graph.graph, first_sequence, aconn, second_sequence, master_graph.velvet_sequences
        )
      gapfilled = true
    end
    scaffold_sequence.gsub!('-','') #remove gaps i.e. where the consensus is a gap
    return scaffold_sequence, variants, gapfilled
  end

  def revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

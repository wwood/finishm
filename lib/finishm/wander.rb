class Bio::FinishM::Wanderer
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :contig_end_length => 200,
    :graph_search_leash_length => 20000,
    :unscaffold_first => false,
    :recoherence_kmer => 1,
    }

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm wander --contigs <contig_file> --fastq-gz <reads..> --output-connections <output.csv>

    Takes a set of contigs/scaffolds from a genome and finds connections in the graph between them. A connection here is given as
    the length of the shortest path between them, without actually computing all the paths.

    This can be used for scaffolding, because if a contig end only connects to one other contig end, then
    those contigs might be scaffolded together.

    This method can also be used for 'pre-scaffolding', in the following sense. If the shortest path between
    two contig ends is 10kb, and a mate pair library with insert size 2kb suggests a linkage
    between the two ends, then the mate pair linkage is likely false (as long as there is sufficient
    coverage in the reads, and not overwhelming amounts of strain heterogeneity, etc.).

    Example:

finishm wander --contigs contigs.fasta --fastq-gz reads.1.fq.gz,reads.2.fq.gz --output-directory finishm_wander_results

    That will create a collapsed de-Bruijn graph from reads.1.fq.gz and reads.2.fq.gz, then try to find connections between
    the starts and the ends of the contigs in contigs.fasta through the de-Bruijn graph. Any connections found
    are output to connections.csv

    \n\n"

    options.merge!(DEFAULT_OPTIONS)

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
      options[:contigs_file] = arg
    end

    optparse_object.separator "\nOutput modes:\n\n"
    optparse_object.on("--output-scaffolds PATH", "Output scaffolds in FASTA format [required]") do |arg|
      options[:output_scaffolds_file] = arg
    end
    optparse_object.on("--output-connections PATH", "Output connections in tab-separated format [required]") do |arg|
      options[:output_connection_file] = arg
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
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs[default: #{options[:unscaffold_first] }]") do
      options[:unscaffold_first] = true
    end
    optparse_object.on("--proceed-on-short-contigs", "By default, when overly short contigs are encountered, finishm croaks. This option stops the croaking [default: #{options[:proceed_on_short_contigs] }]") do
      options[:proceed_on_short_contigs] = true
    end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      [
        :contigs_file,
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end
      if options[:output_scaffolds_file].nil? and
        options[:output_connection_file].nil?
        return "Need to specify either output scaffolds or output connections file"
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv=[])
    # Read in all the contigs sequences, removing those that are too short
    probe_sequences = []
    sequences_to_probe = {}
    overly_short_sequence_count = 0
    process_sequence = lambda do |name, seq|
      if seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{name}"
        overly_short_sequence_count += 1
        nil
      else
        sequences_to_probe[name] = seq

        sequence = seq.seq
        fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
        probe_sequences.push fwd2.reverse_complement.to_s

        probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]

        # 'return' the probe indices that have been assigned
        [probe_sequences.length-2, probe_sequences.length-1]
      end
    end

    scaffolds = nil #Array of Bio::FinishM::ScaffoldBreaker::Scaffold objects. Only set when options[:unscaffold_first] is true
    scaffolded_contig_to_probe_ids = {}
    if options[:unscaffold_first]
      log.info "Unscaffolding scaffolds (before trying to connect them together again)"
      scaffolds = Bio::FinishM::ScaffoldBreaker.new.break_scaffolds options[:contigs_file]
      scaffolds.each do |scaffold|
        scaffold.contigs.each do |contig|
          process_sequence.call contig.name, contig.sequence
        end
      end
    else
      # Else don't split up any of the sequences
      log.info "Reading input sequences.."
      Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
        process_sequence.call seq.definition, seq.seq
      end
    end

    if overly_short_sequence_count > 0
      unless options[:proceed_on_short_contigs]
        raise "Not proceding as some contigs are too short (length < 2 * overhang). You might try: "+
          "(1) omitting the smaller contigs, (2) reducing the --overhang parameter, or "+
          "(3) using --proceed-on-short-contigs to continue optimistically ignoring the #{overly_short_sequence_count} short contigs"
      end
    end

    log.info "Searching from #{probe_sequences.length} different contig ends (#{probe_sequences.length / 2} contigs)"

    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Loop over the ends, trying to make connections from each one
    cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new

    log.info "Finding possible connections with recoherence kmer of #{options[:recoherence_kmer] }"
    first_connections = cartographer.wander(finishm_graph, options[:graph_search_leash_length], options[:recoherence_kmer], finishm_graph.velvet_sequences)
    log.info "Found #{first_connections.length} connections with less distance than the leash length, out of a possible #{probe_sequences.length*(probe_sequences.length-1) / 2}"

    probe_descriptions = []
    sequence_names = sequences_to_probe.keys #Ruby hashes are sorted
    (0...finishm_graph.probe_nodes.length).each do |i|
      desc = Bio::FinishM::ConnectionInterpreter::Probe.new
      if i % 2 == 0
        desc.side = :start
        desc.sequence_name = sequence_names[i / 2]
      else
        desc.side = :end
        desc.sequence_name = sequence_names[(i-1) / 2]
      end
      probe_descriptions.push desc
    end

    # Gather connections ready for output
    distance_calibrator = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    all_connections = []
    first_connections.each do |node_indices, distance|
      calibrated_distance = distance_calibrator.calibrate_distance_accounting_for_probes(
        finishm_graph,
        node_indices[0],
        node_indices[1],
        distance
        )

      # It is possible that a connection just larger than the leash length is returned.
      # weed these out.
      conn = Bio::FinishM::ConnectionInterpreter::Connection.new
      conn.probe1 = probe_descriptions[node_indices[0]]
      conn.probe2= probe_descriptions[node_indices[1]]
      conn.distance = calibrated_distance
      if calibrated_distance > options[:graph_search_leash_length]
        log.debug "Disregarding connection #{conn} because it was ultimately outside the allowable leash length" if log.debug?
      else
        all_connections.push conn
      end
    end

    # Determine scaffolding connections
    interpreter = Bio::FinishM::ConnectionInterpreter.new(all_connections, sequences_to_probe)
    connections = interpreter.doubly_single_contig_connections
    log.info "Found #{connections.length} connections between contigs that can be used for scaffolding" if log.info?

    # Write scaffold file out
    circular_scaffold_names = []
    num_contigs_in_circular_scaffolds = 0
    num_singleton_contigs = 0
    num_scaffolded_contigs = 0
    scaffolds = interpreter.scaffolds(connections)
    scaffolds.each do |scaffold|
      if scaffold.circular?
        circular_scaffold_names.push name
        num_contigs_in_circular_scaffolds += scaffold.contigs.length
      elsif scaffold.contigs.length == 1
        num_singleton_contigs += 1
      else
        num_scaffolded_contigs += scaffold.contigs.length
      end
    end
    log.info "Found #{circular_scaffold_names.length} circular scaffolds encompassing #{num_contigs_in_circular_scaffolds} contigs"
    log.info "#{num_scaffolded_contigs} contigs were incorporated into scaffolds"
    log.info "#{num_singleton_contigs} contigs were not incorporated into any scaffolds"

    unless options[:output_scaffolds_file].nil?
      File.open(options[:output_scaffolds_file],'w') do |scaffold_file|
        scaffolds.each do |scaffold|
          name = scaffold.name
          if scaffold.circular?
            name += ' circular'
          end

          scaffold_file.puts ">#{name}"
          # Output the NA sequence wrapped
          scaffold_file.puts scaffold.sequence.gsub(/(.{80})/,"\\1\n").gsub(/\n$/,'')
        end
      end
    end

    # Write out all connections to the given file if wanted
    unless options[:output_connection_file].nil?
      File.open(options[:output_connection_file], 'w') do |out|
        all_connections.each do |conn|
          out.puts [
            "#{conn.probe1.sequence_name}:#{conn.probe1.side}",
            "#{conn.probe2.sequence_name}:#{conn.probe2.side}",
            conn.distance
            ].join("\t")
        end
      end
    end

    log.info "All done."
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

  class WanderResult
    attr_accessor :connections
  end
end

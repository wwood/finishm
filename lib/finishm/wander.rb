class Bio::FinishM::Wanderer
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :contig_end_length => 200,
    :graph_search_leash_length => 20000,
    :unscaffold_first => false,
    :recoherence_kmer => 1,
    }

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm wander --contigs <contig_file> --fastq-gz <reads..> --output-connections <output.csv> --output-scaffolds <output.fasta>

    Takes a set of contigs/scaffolds from a genome and finds connections in the graph between them. A connection here is given as
    the length of the shortest path between them, without actually computing all the paths.

    This can be used for scaffolding, because if a contig end only connects to one other contig end, then
    those contigs might be scaffolded together.

    This method can also be used for 'pre-scaffolding', in the following sense. If the shortest path between
    two contig ends is 10kb, and a mate pair library with insert size 2kb suggests a linkage
    between the two ends, then the mate pair linkage is likely false (as long as there is sufficient
    coverage in the reads, and not overwhelming amounts of strain heterogeneity, etc.).

    Example:

      finishm wander --contigs contigs.fasta --fastq-gz reads.1.fq.gz,reads.2.fq.gz --output-scaffolds scaffolds.fasta

    That will create a collapsed de-Bruijn graph from reads.1.fq.gz and reads.2.fq.gz, then try to find connections between
    the starts and the ends of the contigs in contigs.fasta through the de-Bruijn graph. The new scaffolds are then
    output to scaffolds.fasta

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
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs. [default: #{options[:unscaffold_first] }]") do
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
    contig_sequences = []
    contig_names = []
    overly_short_sequence_count = 0
    process_sequence = lambda do |name, seq|
      if seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{name}"
        overly_short_sequence_count += 1
        nil
      else
        contig_sequences.push seq.to_s
        contig_names.push name

        sequence = seq.seq
        fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
        probe_sequences.push fwd2.reverse_complement.to_s

        probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]

        # 'return' the probe indices that have been assigned
        [probe_sequences.length-2, probe_sequences.length-1]
      end
    end

    scaffolds = nil #Array of Bio::FinishM::ScaffoldBreaker::Scaffold objects.
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

    log.info "Finding possible connections with recoherence kmer of #{options[:recoherence_kmer] }"
    all_connections = probed_graph_to_connections(finishm_graph, contig_sequences, options)
    log.debug "Finished actual wandering, found #{all_connections.length} connections" if log.debug?

    # Determine scaffolding connections
    interpreter = Bio::FinishM::ConnectionInterpreter.new(all_connections, (0...contig_sequences.length))
    connections = interpreter.doubly_single_contig_connections
    log.debug "Found #{connections.length} connections between contigs that can be used for scaffolding" if log.debug?
    scaffolds = interpreter.scaffolds(connections)

    # Gather some stats
    circular_scaffold_names = []
    num_contigs_in_circular_scaffolds = 0
    num_singleton_contigs = 0
    num_scaffolded_contigs = 0
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
        scaffolds.each_with_index do |scaffold, i|
          name = nil
          if scaffold.contigs.length == 1
            name = "scaffold#{i+1}"
          else
            name = "scaffold#{i+1}"
          end
          if scaffold.circular?
            name += ' circular'
          end

          scaffold_file.puts ">#{name}"
          # Output the NA sequence wrapped
          seq = scaffold.sequence(contig_sequences)
          scaffold_file.puts seq.gsub(/(.{80})/,"\\1\n").gsub(/\n$/,'')
        end
      end
    end

    # Write out all connections to the given file if wanted
    unless options[:output_connection_file].nil?
      File.open(options[:output_connection_file], 'w') do |out|
        all_connections.each do |conn|
          out.puts [
            "#{contig_names[conn.probe1.sequence_index]}:#{conn.probe1.side}",
            "#{contig_names[conn.probe2.sequence_index]}:#{conn.probe2.side}",
            conn.distance
            ].join("\t")
        end
      end
    end

    log.info "All done."
  end

  # Given a probed graph, wander between all the nodes, and then return an
  # instance of Bio::FinishM::ConnectionInterpreter::Scaffold. Required options:
  # * :graph_search_leash_length
  # * :recoherence_kmer
  def probed_graph_to_connections(finishm_graph, contig_sequences, options)
    # Loop over the ends, trying to make connections from each one
    cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new

    first_connections = cartographer.wander(finishm_graph, options[:graph_search_leash_length], options[:recoherence_kmer], finishm_graph.velvet_sequences)
    log.debug "Initially found #{first_connections.length} connections with less distance than the leash length" if log.debug?

    probe_descriptions = []
    (0...finishm_graph.probe_nodes.length).each do |i|
      desc = Bio::FinishM::ConnectionInterpreter::Probe.new
      if i % 2 == 0
        desc.side = :start
        desc.sequence_index = i / 2
      else
        desc.side = :end
        desc.sequence_index = (i-1) / 2
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
      conn.probe2 = probe_descriptions[node_indices[1]]
      conn.distance = calibrated_distance
      if calibrated_distance > options[:graph_search_leash_length]
        log.debug "Disregarding connection #{conn} because it was ultimately outside the allowable leash length" if log.debug?
      else
        all_connections.push conn
      end
    end
    return all_connections
  end
end

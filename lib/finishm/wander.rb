class Bio::FinishM::Wanderer
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm gapfill --contig <contig_file> --fastq-gz <reads..> --output-connections <output.csv>

    Takes a set of reads and a contig that contains gap characters. Then it tries to fill in
    these N characters. It is possible that there is multiple ways to close the gap - in that case
    each is reported. \n\n"

    options.merge!({
      :contig_end_length => 200,
      :graph_search_leash_length => 20000,
      :unscaffold_first => false,
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
      options[:contigs_file] = arg
    end
    optparse_object.on("--output-connections PATH", "Output found paths to this file in fasta format [required]") do |arg|
      options[:output_connection_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: #{options[:graph_search_leash_length]}]") do |arg|
      options[:graph_search_leash_length] = arg
    end
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs[default: #{options[:unscaffold_first]}]") do |arg|
      options[:unscaffold_first] = true
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
        :contigs_file,
        :output_connection_file
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv)
    # Read in all the contigs sequences, removing those that are too short
    probe_sequences = []
    sequence_names = []
    process_sequence = lambda do |name, seq|
      if seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{name}"
        next
      end
      sequence_names.push name

      sequence = seq.seq
      fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
      probe_sequences.push fwd2.reverse_complement.to_s

      probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]
    end

    scaffolds = nil #Array of Bio::FinishM::ScaffoldBreaker::Scaffold objects. Only set when options[:unscaffold_first] is true
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
      Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
        process_sequence.call seq.definition, seq.seq
      end
    end
    log.info "Searching from #{probe_sequences.length} different contig ends (#{probe_sequences.length/2} contigs)"

    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Loop over the ends, trying to make connections from each one
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new

    log.info "Finding possible connections with a depth first search"
    first_connections = cartographer.depth_first_search_with_leash(finishm_graph, options[:graph_search_leash_length])
    log.info "Found #{first_connections.length} connections with less distance than the leash length, out of a possible #{probe_sequences.length*(probe_sequences.length-1)/2}"

    File.open(options[:output_connection_file], 'w') do |out|
      first_connections.each do |node_indices, distance|
        sequence_names_and_directions = node_indices.collect do |i|
          node = finishm_graph.probe_nodes[i]
          node_forward = finishm_graph.probe_node_directions[i]

          sequence_name = nil
          side = nil
          if i % 2 == 0
            side = 'start'
            sequence_name = sequence_names[i/2]
          else
            side = 'end'
            sequence_name = sequence_names[(i-1)/2]
          end
          [sequence_name, side]
        end
        out.puts [
          sequence_names_and_directions,
          distance
        ].flatten.join("\t")
      end
    end

    # While this may not happen much, when scaffolding a genome it would be good to know if there are any nodes that:
    # 1. have no connections - possibly this indicates missing assembly
    # 2. have exactly 1 connection - this probably indicates that the connection is the 'true' in the genome
    # 3. if connections from #2 are true, then this may show that other connections are false, because
    #    each node can only be connected to at most 1 other node, and this is true for both ends of the true edge.
    #    After removing those true edges

    #TODO: implemented this in repo hamiltonian_cycler, need to incorporate it here. See also a script in luca/bbbin that uses that library.
    #TODO: look for hamiltonian paths as well as hamiltonian cycles

  end
end

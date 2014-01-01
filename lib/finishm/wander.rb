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
    Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
      if seq.seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{seq.definition}"
        next
      end
      sequence_names.push seq.definition

      sequence = seq.seq
      fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
      probe_sequences.push fwd2.reverse_complement.to_s

      probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]
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

    first_connections.each do |node_indices, distance|
      sequence_names_and_directions = node_indices.collect do |i|
        node = finishm_graph.probe_nodes[i]
        node_forward = finishm_graph.probe_node_directions[i]

        sequence_name = nil
        side = nil
        if i % 2 == 0
          side = 'start'
          sequence_name = sequence_names[start_probe_index/2]
        else
          side = 'end'
          sequence_name = sequence_names[(start_probe_index-1)/2]
        end
        [sequence_name, side]
      end
      puts [
        sequence_names_and_directions,
        distance
      ].flatten.join("\t")
    end
  end
end

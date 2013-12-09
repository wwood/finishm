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
    num_total_connections = 0
    probe_sequences.each_with_index do |probe, start_probe_index|
      start_node = finishm_graph.probe_nodes[start_probe_index]
      start_node_forward = finishm_graph.probe_node_directions[start_probe_index]

      start_node_sequence_name = nil
      start_node_side = nil
      if start_probe_index % 2 == 0
        start_node_side = 'start'
        start_node_sequence_name = sequence_names[start_probe_index/2]
      else
        start_node_side = 'end'
        start_node_sequence_name = sequence_names[(start_probe_index-1)/2]
      end

      if start_node
        probe_sequences.each_with_index do |probe2, terminal_probe_index|
          # Don't try to connect an end with itself
          break if start_probe_index == terminal_probe_index

          end_node = finishm_graph.probe_nodes[terminal_probe_index]
          if end_node
            paths = cartographer.find_trails_between_nodes(
              finishm_graph.graph,
              start_node,
              end_node,
              options[:graph_search_leash_length],
              start_node_forward
            )

            terminal_node_sequence_name = nil
            terminal_node_side = nil
            if terminal_probe_index % 2 == 0
              terminal_node_side = 'start'
              terminal_node_sequence_name = sequence_names[terminal_probe_index/2]
            else
              terminal_node_side = 'end'
              terminal_node_sequence_name = sequence_names[(terminal_probe_index-1)/2]
            end
            log.debug "Trying to connect #{start_node_sequence_name}/#{start_node_side} to #{terminal_node_sequence_name}/#{terminal_node_side}" if log.debug?

            unless paths.empty?
              num_total_connections += 1

              puts [
                start_node_sequence_name,
                start_node_side,
                terminal_node_sequence_name,
                terminal_node_side,
                paths.length,
              ].join("\t")
            end
          else
            log.warn "Unable to retrieve probe sequences for probe index #{start_probe_index}, so cannot find connections towards this contig side"
          end
        end
      else
        #TODO: the below error message should link up better with the contig names, to help users a bit more. Or maybe it is unlikely to happen much any more anyway.
        log.warn "Unable to retrieve probe sequences for probe index #{start_probe_index}, so cannot find connections from this contig side"
      end
    end
    log.info "Found #{num_total_connections} gap filling(s) in total, out of a possible #{probe_sequences.length*(probe_sequences.length-1)/2}"
  end
end

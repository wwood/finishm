class Bio::FinishM::Fluff
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm fluff --contigs <contig_file> --fastq-gz <reads..> --output-fluff-file <output.csv>

    Takes a set of contigs, and places probes across them (e.g. every 2kb), and then explores the
    graph from each of these probes, taking all paths within some leash length, including the 'fluff'
    which is not the same path as along the contig. Prints out all of these paths to a fasta file.\n\n"

    options.merge!({
      :probe_spacing => 2000,
      :probe_length => 100,
      :graph_search_leash_length => 20000,
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "fasta file containing contigs to find the fluff on [required]") do |arg|
      options[:contigs_file] = arg
    end
    optparse_object.on("--output-fluff-file PATH", "Output found paths to this file in fasta format [required]") do |arg|
      options[:output_fluff_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--probe-spacing NUM", Integer, "Distance between probe points in the contig [default: #{options[:probe_spacing]}]") do |arg|
      options[:probe_spacing] = arg
    end
    optparse_object.on("--probe-size NUM", Integer, "Length of the probe to be inserted into the velvet graph. Must be greater than graph kmer length. [default: #{options[:probe_length]}]") do |arg|
      options[:probe_length] = arg
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
        :output_fluff_file
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      unless options[:velvet_kmer_size] < options[:probe_length]
        return "The probe length must be greater than the kmer length, otherwise it will not be incorporated into the kmer graph"
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv)
    # Read in all the contigs sequences
    probe_sequences = []
    sequence_names = []
    Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
      sequence_names.push seq.definition

      sequence = seq.seq
      0.step(sequence.length-1-options[:probe_length], options[:probe_spacing]) do |offset|
        # Only probe in the forward direction
        probe_sequence = sequence[offset...offset+options[:probe_length]]
        probe_sequences.push probe_sequence
      end
    end
    log.info "Searching from #{probe_sequences.length} different probes from (#{sequence_names.length}) contigs)"

    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Loop over the ends, trying to make connections from each one
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    fluffings = fluffer.fluff(finishm_graph)

    log.info "Found #{num_total_connections} gap filling(s) in total, out of a possible #{probe_sequences.length*(probe_sequences.length-1)/2}"

    # Print out the sequences
    File.open(options[:output_fluff_file], 'w') do |output|
      fluffings.each_with_index do |path_set, probe_number|
        path_set.each do |path, path_number|
          output.puts ">probe#{probe_number+1}_path#{path_number+1}"
          output.puts path.sequence
        end
      end
    end
  end
end

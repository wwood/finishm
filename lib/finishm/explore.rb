class Bio::FinishM::Explorer
  include Bio::FinishM::Logging

  class InterestingPlace
    attr_accessor :contig_name, :start_or_end

    def to_s
      "#{@contig_name}:#{@start_or_end}"
    end
  end

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm explore --contigs <contig_file> --interesting-ends <contig1:end,contig2:start,..> --fastq-gz <reads..> --output-explored-paths <output.fa>

    Given a contig end, explore the assembly graph to determine what is out there. Does assembly
    fail because of lack of coverage? Is there more sequence out there that has yet to be explored?
    \n\n"

    options.merge!({
      :contig_end_length => 200,
      :graph_search_leash_length => 20000,
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "Fasta file containing contigs to find the fluff on [required]") do |arg|
      options[:contigs_file] = arg
    end
    optparse_object.on("--interesting-ends INTERESTING_PLACES", Array, "Comma-separated list of places to explore from e.g. 'contig1:end,MyContig2:start' to explore from the end of contig1 and the start of MyContig2. Names of contigs are as they are in the given --contigs file. [required]") do |arg|
      arg.each do |tuple|
        options[:interesting_places] ||= []
        splits = tuple.split(':')
        if splits.length != 2
          log.error "Unable to parse this --interesting-ends argument: #{tuple}"
          exit 1
        end
        place = InterestingPlace.new
        place.contig_name = splits[0]
        if %(start end).include?(splits[1])
          place.start_or_end = splits[1]
        else
          log.error "Unable to parse this --interesting-ends argument, second half must be 'start' or 'end': #{tuple}"
          exit 1
        end
        options[:interesting_places].push place
      end
    end
    optparse_object.on("--output-explored-paths PATH", "Output found paths to this file in fasta format [required]") do |arg|
      options[:output_trails_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs[default: #{options[:graph_search_leash_length]}]") do |arg|
      options[:unscaffold_first] = true
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
        :interesting_places,
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
        if options[:contig_end_length] < options[:velvet_kmer_size]
          return "The overhang must be greater than the size of the assembly kmer"
        else
          return nil
        end
      else
        return nil
      end
    end
  end

  def run(options, argv)
    # Read in all the contigs sequences, removing those that are too short
    probe_sequences = []
    sequence_names = []
    interesting_place_probe_indices = []
    process_sequence = lambda do |name,seq|
      if seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{name}"
        next
      end
      if sequence_names.include?(name)
        log.error "Found duplicate sequence names, being conservative and not continuuing #{name}"
        exit 1
      end
      sequence_names.push name

      sequence = seq.seq
      fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
      probe_sequences.push fwd2.reverse_complement.to_s

      probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]
    end

    scaffolds = nil
    if options[:unscaffold_first]
      log.info "Unscaffolding scaffolds (before trying to connect them together again)"
      scaffolds = Bio::FinishM::ScaffoldBreaker.new.break_scaffolds options[:contigs_file]
      scaffolds.each do |scaffold|
        scaffold.contigs.each do |contig|
          process_sequence.call contig.name, contig.sequence
        end
      end
    else
      Bio::FlatFile.foreach(options[:contigs_file]) do |s|
        process_sequence.call s.definition, s.seq
      end
    end

    # Collect the node IDs that I'm interested in before generating the graph so don't have to do a whole assembly before getting an argument error
    interesting_probe_ids_to_place = {}
    options[:interesting_places].each do |place|
      seq_index = sequence_names.find_index place.contig_name
      if seq_index.nil?
        log.error "Unable to find interesting contig #{place.contig_name}, cannot continue"
        exit 1
      else
        base = seq_index*2
        if place.start_or_end == 'start'
          #
        elsif place.start_or_end == 'end'
          base += 1
        end
        interesting_probe_ids_to_place[base] = place
      end
    end


    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Explore from the interesting nodes

    output = nil
    if options[:output_trails_file] == '-'
      log.info "When trails are found, writing them to stdout"
      output = $stdout
    else
      log.info "When trails are found, writing them to #{options[:output_trails_file]}"
      output = File.open(options[:output_trails_file],'w')
    end

    explorer = Bio::AssemblyGraphAlgorithms::GraphExplorer.new
    interesting_probe_ids_to_place.each do |probe_id, place|
      log.info "Exploring from #{place}"
      if finishm_graph.probe_nodes[probe_id].nil?
        log.warn "Unable to find anchor node for #{place}, skipping exploration from there"
      else
        # Do exploration
        onode = finishm_graph.initial_path_from_probe probe_id
        paths = explorer.explore_from_node(finishm_graph.graph, onode, options[:graph_search_leash_length])
        log.info "Found #{paths.length} paths from #{place}"

        # Print explorations that come back
        paths.each do |explore_path|
          output.puts ">#{place.contig_name}:#{place.start_or_end} #{explore_path.termination_type} nodes:#{explore_path}"
          output.puts explore_path.path.sequence
        end
      end
    end
    output.close unless output == $stdout
  end
end

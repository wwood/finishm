class Bio::FinishM::GapFiller
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm gapfill --contigs <contigs_file> --fastq-gz <reads..> --output-trails-fasta <output.fa>

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
    optparse_object.on("--output-trails-fasta PATH", "Output overview of connections between contigs to this file [required]") do |arg|
      options[:overall_trail_output_fasta_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)
    optparse_object.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
      options[:output_graph_png] = arg
    end
    optparse_object.on("--assembly-svg PATH", "Output assembly as an SVG file [default: off]") do |arg|
      options[:output_graph_svg] = arg
    end
    optparse_object.on("--assembly-dot PATH", "Output assembly as an DOT file [default: off]") do |arg|
      options[:output_graph_dot] = arg
    end

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", "Start assembling this far from the gap [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
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
        :overall_trail_output_fasta_file
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}"
        end
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv)
    # Read in all the contigs sequences and work out where the gaps are
    sequences = {}
    gaps = []
    #TODO: refactor this code so that it uses ScaffoldBreaker, because DRY
    Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
      seq_name = seq.definition
      sequences[seq_name] = seq.seq

      unless seq.seq.match(/^[ATGCN]+$/i)
        log.warn "Found unexpected characters in the sequence, continuing optimistically, but not quite sure what will happen.. good luck"
      end

      # Find all Ns in the current sequence
      seq.seq.scan(/N+/i) do
        gap = Gap.new
        gap.parent_sequence_name = seq_name
        gap.start = $~.offset(0)[0]
        gap.stop = $~.offset(0)[1]
        if gap.start == 0 or gap.stop == sequences[seq_name].length
          log.warn "Ignoring gap #{gap.coords} because it is too close to one (or both) end(s) of the contig"
        else
          gaps.push gap
        end
      end
    end
    log.info "Detected #{gaps.length} gaps from #{gaps.collect{|g| g.parent_sequence_name}.uniq.length} different sequence(s)"

    probe_sequences = gaps.collect do |gap|
      sequence = sequences[gap.parent_sequence_name]

      if gap.start < options[:contig_end_length] or gap.stop > sequence.length - options[:contig_end_length]
        log.warn "Found a gap that was too close to the end of a contig - you might have problems: #{gap.coords}"
      end

      fwd2 = Bio::Sequence::NA.new(sequence[gap.stop..(gap.stop+options[:contig_end_length])])
      probes = [
        sequence[gap.start-options[:contig_end_length]-1...gap.start],
        fwd2.reverse_complement.to_s
      ]
      #TODO: this could probably be handled better.. e.g. make sure the kmer is less than the amount of seq on the end, etc.
      if probes[0].match(/N/i) or probes[1].match(/N/i)
        log.warn "Noticed gap that was too close together, not sure if things will work out for #{gap.coords}"
      end
      probes
    end.flatten
    log.debug "Generated #{probe_sequences.length} probes e.g. #{probe_sequences[0]}"

    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Output optional graphics. TODO: these should be abstracted because it is not DRY
    if options[:output_graph_png] or options[:output_graph_svg] or options[:output_graph_dot]
      viser = Bio::Assembly::ABVisualiser.new
      # TODO: make these visualise more than one join somehow
      gv = viser.graphviz(finishm_graph.graph, {
        :start_node_id => finishm_graph.probe_nodes[0].node_id,
        :end_node_id => finishm_graph.probe_nodes[1].node_id})

      if options[:output_graph_png]
        log.info "Converting assembly to a graphviz PNG"
        gv.output :png => options[:output_graph_png], :use => :neato
      end
      if options[:output_graph_svg]
        log.info "Converting assembly to a graphviz SVG"
        gv.output :svg => options[:output_graph_svg], :use => :neato
      end
      if options[:output_graph_dot]
        log.info "Converting assembly to a graphviz DOT"
        gv.output :dot => options[:output_graph_dot]
      end
    end

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    num_total_trails = 0
    File.open(options[:overall_trail_output_fasta_file],'w') do |f|
      log.info "Searching for trails between the nodes within the assembly graph"
      (0...(probe_sequences.length / 2)).collect{|i| i*2}.each do |start_probe_index|
        gap_number = start_probe_index / 2
        gap = gaps[gap_number]
        log.info "Now working through gap number #{gap_number}: #{gap.coords}"
        start_node = finishm_graph.probe_nodes[start_probe_index]
        start_node_forward = finishm_graph.probe_node_directions[start_probe_index]
        end_node = finishm_graph.probe_nodes[start_probe_index+1]
        if start_node and end_node
          trails = cartographer.find_trails_between_nodes(finishm_graph.graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
          log.info "Found #{trails.length} trails for #{gap.coords}"
          trails.each_with_index do |trail, i|
            #TODO: need to output this as something more sensible e.g. VCF format
            f.puts ">#{gap.coords}_trail#{i+1}"
            f.puts trail.sequence
            num_total_trails += 1
          end
          # TODO: If there is only 1 trail, then output scaffolding information

        else
          log.warn "Unable to retrieve both probes from the graph for gap #{gap_number} (#{gap.coords}), skipping"
        end
      end
    end
    log.info "Found #{num_total_trails} gap filling(s) in total"
  end

  class Gap
    attr_accessor :parent_sequence_name, :start, :stop

    def coords
      @parent_sequence_name.to_s+':'+(@start+1).to_s+'-'+(@stop).to_s
    end
  end
end

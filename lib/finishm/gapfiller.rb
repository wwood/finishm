class Bio::FinishM::GapFiller
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm gapfill --contigs <contigs_file> --fastq-gz <reads..> --output-trails-fasta <output.fa>

Takes a set of reads and a contig that contains gap characters. Then it tries to fill in
these N characters. It is possible that there is multiple ways to close the gap - in that case
each can be reported. \n\n"

    options.merge!({
      :contig_end_length => 200,
      :graph_search_leash_length => 20000,
      })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
      options[:contigs_file] = arg
    end
    optparse_object.on("--output-fasta PATH", "Output the gap-filled sequence to this file [required]") do |arg|
      options[:overall_fasta_file] = arg
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
    optparse_object.on("--overhang NUM", "Start assembling this far from the gap [default: #{options[:contig_end_length] }]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--output-trails-fasta PATH", "Output all connections between contigs to this file [default: don't output]") do |arg|
      options[:overall_trail_output_fasta_file] = arg
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
        :overall_fasta_file
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
    scaffolds = Bio::FinishM::ScaffoldBreaker.new.break_scaffolds(options[:contigs_file])
    gaps = []
    output_fasta_file = File.open(options[:overall_fasta_file],'w')
    num_without_gaps = 0
    scaffolds.each do |scaffold|
      sgaps = scaffold.gaps
      if sgaps.empty?
        num_without_gaps += 1
        output_fasta_file.puts ">#{scaffold.name }"
        output_fasta_file.puts scaffold.sequence
      else
        gaps.push scaffold.gaps
      end
    end
    gaps.flatten!
    log.info "Detected #{gaps.length} gaps from #{scaffolds.length} different sequence(s). #{num_without_gaps } sequences were gap-free"

    # Create probe sequences
    probe_sequences = gaps.collect do |gap|
      sequence = gap.scaffold.sequence

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
    log.debug "Generated #{probe_sequences.length} probes e.g. #{probe_sequences[0] }"


    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Output optional graphics.
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


    # Do the gap-filling and print out the results
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
    num_total_trails = 0
    num_singly_filled = 0
    num_unbridgable = 0

    output_trails_file = nil
    output_trails_file = File.open(options[:overall_trail_output_fasta_file],'w') unless options[:overall_trail_output_fasta_file].nil?

    # Print the fasta output for the scaffold
    print_scaffold = lambda do |last_scaffold, gapfilled_sequence|
      output_fasta_file.puts ">#{last_scaffold.name }"
      gapfilled_sequence.push last_scaffold.contigs[last_scaffold.contigs.length-1] #add last contig
      output_fasta_file.puts gapfilled_sequence
    end
    # Lambda to add a gap the the String representing the scaffold
    filler = lambda do |trails, following_contig|
      if trails.length == 1
        # If there is only 1 trail, then output scaffolding information
        acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
        acon.start_probe_node = start_node
        acon.end_probe_node = end_node
        acon.start_probe_read_id = start_probe_index+1
        acon.end_probe_read_id = start_probe_index+2
        acon.start_probe_contig_offset = options[:contig_end_length]
        acon.end_probe_contig_offset = options[:contig_end_length]
        acon.paths = trails
        gapfilled_sequence = printer.one_connection_between_two_contigs(
          finishm_graph.graph,
          gapfilled_sequence,
          following_contig.sequence
          )
        num_singly_filled += 1
      else
        # Otherwise don't make any assumptions
        num_unbridgable += 1 if trails.empty?
        # TODO: even the there is multiple trails, better info can still be output here
        gapfilled_sequence += 'N'*gap.length
      end
    end

    log.info "Searching for trails between the nodes within the assembly graph"
    gapfilled_sequence =
    last_scaffold = nil
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

        # print the sequences of the trails if asked for:
        trails.each_with_index do |trail, i|
          #TODO: need to output this as something more sensible e.g. VCF format
          unless output_trails_file.nil?
            output_trails_file.puts ">#{gap.coords}_trail#{i+1}"
            output_trails_file.puts trail.sequence
          end
          num_total_trails += 1
        end

        # Output the updated sequence. Fill in the sequence if there is only 1 trail
        if gap.scaffold == last_scaffold
          gapfilled_sequence.push scaffold.contigs[gap.number]
          filler.call trails, scaffold.contigs[gap.number+1]
        else
          unless last_scaffold.nil?
            # print the gapfilled (or not) scaffold.
            print_scaffold.call(last_scaffold, gapfilled_sequence)
          end
          #reset
          last_scaffold = scaffold

          #add the current gap (and the contig before it)
          gapfilled_sequence = scaffold.contigs[gap.number]
          filler.call trails, scaffold.contigs[gap.number+1]
        end
      else
        raise "Unable to retrieve both probes from the graph for gap #{gap_number} (#{gap.coords}), fail"
      end
    end
    print_scaffold.call(last_scaffold, gapfilled_sequence) # print the last scaffold

    log.info "#{num_unbridgable } had no bridging paths in the graph within the leash, and found #{num_total_trails} trails in total."
    log.info "Filled #{num_singly_filled } out of #{gaps.length } gaps."

    output_trails_file.close unless output_trails_file.nil?
    output_fasta_file.close
  end
end

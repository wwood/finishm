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

    optparse_object.separator "\nVisualisation options (of all joins):\n\n" #TODO improve this help
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
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the gap [default: #{options[:contig_end_length] }]") do |arg|
      options[:contig_end_length] = arg
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this many base pairs and not much more [default: #{options[:graph_search_leash_length] }]") do |arg|
      options[:graph_search_leash_length] = arg
    end
    optparse_object.on("--disambiguation-kmer NUM", Integer, "Use a kmer longer than the original velvet one, to help remove bubbles and circular paths [default: none] }]") do |arg|
      options[:recoherence_kmer] = arg
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
    log.info "Detected #{gaps.length} gap(s) from #{scaffolds.length} different contig(s). #{num_without_gaps } contig(s) were gap-free."

    # Create probe sequences
    probe_sequences = gaps.collect do |gap|
      sequence = gap.scaffold.sequence

      if gap.start < options[:contig_end_length] or gap.stop > sequence.length - options[:contig_end_length]
        log.warn "Found a gap that was too close to the end of a contig - you might have problems: #{gap.coords}"
      end

      log.debug "Processing gap number #{gap.number}, #{gap.coords}"
      first_coords = [
        gap.start-options[:contig_end_length]-1,
        gap.start-1,
        ]
      second_coords = [
        gap.stop,
        (gap.stop+options[:contig_end_length]),
        ]
      log.debug "Coordinates of the probes are #{first_coords} and #{second_coords}"
      second = sequence[second_coords[0]..second_coords[1]]
      probes = [
        sequence[first_coords[0]...first_coords[1]],
        Bio::Sequence::NA.new(second).reverse_complement.to_s,
        ]
      #TODO: this could probably be handled better.. e.g. if the amount of sequence is too small, just throw it out and make one big gap
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
      #gapfilled_sequence += last_scaffold.contigs[last_scaffold.contigs.length-1].sequence #add last contig
      output_fasta_file.puts gapfilled_sequence
    end
    # Lambda to add a gap the the String representing the scaffold
    #TODO: if the trail is not filled then the wrong sequence is currently printed. BUG???
    filler = lambda do |trails, following_contig, start_node, end_node_inward, start_probe_index, gapfilled_sequence, gap|
      gapfilled = nil
      if trails.length == 1
        # If there is only 1 trail, then output scaffolding information
        acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
        acon.start_probe_node = start_node.node
        acon.end_probe_node = end_node_inward.node
        acon.start_probe_read_id = start_probe_index+1
        acon.end_probe_read_id = start_probe_index+2
        acon.start_probe_contig_offset = options[:contig_end_length]
        acon.end_probe_contig_offset = options[:contig_end_length]
        acon.paths = trails
        gapfilled = printer.one_connection_between_two_contigs(
          finishm_graph.graph,
          gapfilled_sequence,
          acon,
          following_contig.sequence
          )
        num_singly_filled += 1
        #puts ">gapfilled"
        #puts gapfilled_sequence
      else
        # Otherwise don't make any assumptions
        num_unbridgable += 1 if trails.empty?
        # TODO: even the there is multiple trails, better info can still be output here
        gapfilled = gapfilled_sequence + 'N'*gap.length + following_contig.sequence
      end
      gapfilled #return this string
    end

    log.info "Searching for trails between the nodes within the assembly graph"
    log.info "Using contig overhang length #{options[:contig_end_length]} and leash length #{options[:graph_search_leash_length]}"
    gapfilled_sequence = ''
    last_scaffold = nil
    (0...(probe_sequences.length / 2)).collect{|i| i*2}.each do |start_probe_index|
      gap_number = start_probe_index / 2
      gap = gaps[gap_number]
      log.info "Now working through gap number #{gap_number}: #{gap.coords}"
      start_onode = finishm_graph.velvet_oriented_node(start_probe_index)
      end_onode_inward = finishm_graph.velvet_oriented_node(start_probe_index+1)
      if start_onode and end_onode_inward
        # The probe from finishm_graph points in the wrong direction for path finding
        end_onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
        end_onode.node = end_onode_inward.node
        end_onode.first_side = end_onode_inward.starts_at_start? ? Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST : Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

        trails = cartographer.find_trails_between_nodes(
          finishm_graph.graph, start_onode, end_onode, options[:graph_search_leash_length],
          :recoherence_kmer => options[:recoherence_kmer]
        )
        log.info "Found #{trails.trails.length} trails for #{gap.coords}"
        if trails.circular_paths_detected
          log.warn "Circular path detected here, not attempting to gapfill"
        end
        # Convert the trails into OrientedNodePaths
        trails = trails.collect do |trail|
          path = Bio::Velvet::Graph::OrientedNodeTrail.new
          path.trail = trail
          path
        end

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
          # We are still building the current scaffold
          #gapfilled_sequence += gap.scaffold.contigs[gap.number].sequence
          log.debug "Before adding next chunk of contig, length of scaffold being built is #{gapfilled_sequence.length}" if log.debug?
          gapfilled_sequence = filler.call trails, gap.scaffold.contigs[gap.number+1], start_onode, end_onode_inward, start_probe_index, gapfilled_sequence, gap
          log.debug "After adding next chunk of contig, length of scaffold being built is #{gapfilled_sequence.length}" if log.debug?
        else
          # We are onto a new scaffold. Print the previous one (unless this the first one)
          unless last_scaffold.nil?
            # print the gapfilled (or not) scaffold.
            print_scaffold.call(last_scaffold, gapfilled_sequence)
          end
          #reset
          last_scaffold = gap.scaffold

          #add the current gap (and the contig before it)
          log.debug "Before adding first chunk of contig, length of scaffold being built is #{gapfilled_sequence.length}"
          gapfilled_sequence = gap.scaffold.contigs[gap.number].sequence
          log.debug "After adding first chunk of contig, length of scaffold being built is #{gapfilled_sequence.length}"
          gapfilled_sequence = filler.call trails, gap.scaffold.contigs[gap.number+1], start_onode, end_onode_inward, start_probe_index, gapfilled_sequence, gap
          log.debug "After adding first gap sequence and next contig, gapfilled sequence length is #{gapfilled_sequence.length}"
        end
      else
        raise "Unable to retrieve both probes from the graph for gap #{gap_number} (#{gap.coords}), fail"
      end
    end
    print_scaffold.call(last_scaffold, gapfilled_sequence) # print the last scaffold

    log.info "#{num_unbridgable } gaps had no suitable bridging paths in the graph within the leash, and found #{num_total_trails} trails in total."
    log.info "Filled #{num_singly_filled } out of #{gaps.length } gaps."

    output_trails_file.close unless output_trails_file.nil?
    output_fasta_file.close
  end
end

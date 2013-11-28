class Bio::FinishM::Finisher
  include Bio::FinishM::Logging

  def add_options(opts, options)
    options.merge!({
      :min_leftover_length => false,
      :kmer_coverage_target => 1,
      :contig_end_length => 300,
      :graph_search_leash_length => 20000,
      :reads_to_assemble => nil,
    })
    opts.banner = "finishm finish <options>\n\n"

    # TODO: make a better interface for this. Maybe specify an entire genome, and then "Contig_1 end, Contig_3 start" or something
    # Look at the last 300bp of the first contig.
    extract_exactly_one_contig_from_file = lambda do |fasta_file_path|
      contig = nil
      Bio::FlatFile.foreach(Bio::FastaFormat, fasta_file_path) do |e|
        if contig.nil?
          contig = e.seq
        else
          raise "Multiple sequences found in a contig file! I need exactly one"
        end
      end
      raise "I need a contig to be in the start contig file" if contig.nil?
      Bio::Sequence::NA.new(contig.to_s)
    end

    opts.on("--pattern PATTERN", "kmer abundance pattern e.g. '0111001110' [required]") do |arg|
      options[:pattern] = arg
    end
    opts.on("--kmer-abundances FILE", "kmer multiple abundance file [required]") do |arg|
      options[:kmer_multiple_abundance_file] = arg
    end
    opts.on("--upper-threshold NUM", "kmer frequency cutoff to saying 'present' [required]") do |arg|
      options[:upper_threshold] = arg.to_i
    end
    opts.on("--lower-threshold NUM", "kmer frequency cutoff to saying 'not present' [required]") do |arg|
      options[:lower_threshold] = arg.to_i
    end
    opts.on("--reads FILES", "comma-separated list of sequence reads files in the same order as the pattern was supplied [required]") do |arg|
      options[:reads_files] = arg.split(',').collect{|r| File.absolute_path r}
    end
    opts.on("--start-contig FASTA", "path to a fasta file with the starting contig in it (only). Assumes we are building off the end of this contig [required]") do |arg|
      options[:start_contig] = extract_exactly_one_contig_from_file.call arg
    end
    opts.on("--end-contig FASTA", "path to a fasta file with the ending contig in it (only). Assumes we are building onto the start of this contig [required]") do |arg|
      options[:end_contig] = extract_exactly_one_contig_from_file.call arg
    end

    opts.separator "\nOptional arguments:\n\n"
    opts.on("--min-leftover-read-length NUMBER", "when searching for reads with kmers, require the kmer to be at the beginning or end of the selected read [default: #{options[:min_leftover_length]}]") do |arg|
      options[:min_leftover_length] = arg.to_i
    end
    opts.on("--kmer-coverage-target NUMBER", "when searching for reads with kmers, require this many copies per kmer [default: #{options[:kmer_coverage_target]}]") do |arg|
      options[:kmer_coverage_target] = arg.to_i
    end
    opts.on("--already-patterned-reads FILE", "Attempt to assemble the reads in the specified file, useful for re-assembly [default: off]") do |arg|
      options[:already_patterned_reads] = arg
    end
    opts.on("--output-assembly PATH", "Output assembly intermediate files to this directory [default: off]") do |arg|
      options[:output_assembly_path] = arg
    end
    opts.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
      options[:output_graph_png] = arg
    end
    opts.on("--assembly-svg PATH", "Output assembly as an SVG file [default: off]") do |arg|
      options[:output_graph_svg] = arg
    end
    opts.on("--assembly-dot PATH", "Output assembly as an DOT file [default: off]") do |arg|
      options[:output_graph_dot] = arg
    end
    opts.on("--assembly-coverage-cutoff NUMBER", "Require this much coverage in each node, all other nodes are removed [default: #{options[:assembly_coverage_cutoff]}]") do |arg|
      options[:assembly_coverage_cutoff] = arg.to_f
    end
    opts.on("--contig-end-length LENGTH", "Number of base pairs to start into the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
    end

    Bio::FinishM::GraphGenerator.new.add_options opts, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      [:upper_threshold,
      :lower_threshold,
      :pattern,
      :kmer_multiple_abundance_file,
      :reads_files].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}"
        end
      end
    end
    return nil #if here, options all were parsed successfully
  end

  #TODO: this method is too long - split it up by refactoring
  def run(options, argv)
    pooled_reads_filename = 'pooled_sampled_reads.fasta'
    if options[:already_patterned_reads] #If skipping read extraction
      pooled_reads_filename = options[:already_patterned_reads]

    else
      # Parse pattern from cmdline
      desired_pattern = KmerAbundancePattern.new
      desired_pattern.parse_from_human(options[:pattern])
      if options[:reads_files].length != desired_pattern.length
        raise "Number of entries in the pattern #{desired_pattern.length} and number of reads files #{options[:reads].length} not equivalent!"
      end

      # Collect the kmers that will be used to find trusted reads i.e.
      # Go through each line of the kmer abundance file, looking for kmers that suit the pattern
      input_file = File.open options[:kmer_multiple_abundance_file]
      csv = CSV.new(input_file, :col_sep => ' ')

      whitelist_kmers = []
      blacklist_kmers = []
      csv.each do |row|
        max_i = row.length - 2 if max_i.nil?

        kmer = row[0]
        counts = row[1...row.length].collect{|s| s.to_i}
        this_pattern = []
        counts.each_with_index do |count, i|
          if count > options[:upper_threshold]
            this_pattern[i] = true
          elsif count < options[:lower_threshold]
            this_pattern[i] = false
          else
            # coverage was in no man's land between thresholds.
            # Ignore this kmer as noise.
            this_pattern[i] = '-'
          end
        end
        #log.debug "Found pattern #{this_pattern} from kmer #{kmer}, which has abundances #{counts}" if log.debug?

        if desired_pattern.consistent_with? this_pattern
          whitelist_kmers.push row[0]
        else
          # kmer is not present when it should be
          blacklist_kmers.push row[0]
        end
      end
      log.info "After parsing the kmer multiple abundance file, found #{whitelist_kmers.length} kmers that matched the pattern, and #{blacklist_kmers.length} that didn't"
      unless whitelist_kmers.length > 0
        log.error "No kmers found that satisfy the given pattern, exiting.."
        exit 1
      end


      #outdir = options[:output_directory]
      #Dir.mkdir outdir unless Dir.exist?(outdir)

      # grep the pattern out from the raw reads, subsampling so as to not overwhelm the assembler
      #Tempfile.open('whitelist') do |white|
      File.open 'whitelist', 'w' do |white|
        white.puts whitelist_kmers.join("\n")
        white.close

        #Tempfile.open('blacklist') do |black|
        File.open('black','w') do |black|
          black.puts blacklist_kmers.join("\n")
          black.close

          threadpool = []
          sampled_read_files = []
          log.info "Extracting reads that contain suitable kmers"
          options[:reads_files].each_with_index do |file, i|
            next unless desired_pattern[i] #Don't extract reads from reads where those reads should not have been amplified

            sampled = File.basename(file)+'.sampled_reads.fasta'
            sampled_read_files.push sampled

            grep_path = "#{ENV['HOME']}/git/priner/bin/read_selection_by_kmer "
            if options[:min_leftover_length]
              grep_path += "--min-leftover-length #{options[:min_leftover_length]} "
            end
            thr = Thread.new do
              grep_cmd = "#{grep_path} --whitelist #{white.path} --blacklist #{black.path} --reads #{file} --kmer-coverage-target #{options[:kmer_coverage_target]} > #{sampled}"
              log.debug "Running cmd: #{grep_cmd}"
              status, stdout, stderr = systemu grep_cmd
              log.debug stderr

              raise unless status.exitstatus == 0
              log.debug "Finished extracting reads from #{file}"
            end
            threadpool.push thr
          end
          threadpool.each do |thread| thread.join; end #wait until everything is finito

          log.info "Finished extracting reads for sampling. Now pooling sampled reads"
          pool_cmd = "cat #{sampled_read_files.join ' '} >#{pooled_reads_filename}"
          log.debug "Running cmd: #{pool_cmd}"
          status, stdout, stderr = systemu pool_cmd
          raise stderr if stderr != ''
          raise unless status.exitstatus == 0
        end
      end
    end

    log.info "Extracting dummy reads from the ends of contigs to use as anchors"
    start_contig = options[:start_contig]
    end_contig = options[:end_contig]
    if [start_contig.length, end_contig.length].min < 2*options[:contig_end_length]
      log.warn "Choice of initial/terminal nodes to perform graph search with may not be optimal due to the small contig size"
    end
    if [start_contig.length, end_contig.length].min < options[:contig_end_length]
      log.error "At least one contig too small to proceed with current code base, need to fix the code to allow such a small contig"
      exit 1
    end

    probe_sequences = [
      start_contig[start_contig.length-options[:contig_end_length]...start_contig.length],
      Bio::Sequence::NA.new(end_contig[0...options[:contig_end_length]]).reverse_complement.to_s
    ]
    read_input = Bio::FinishM::ReadInput.new
    read_input.fasta_singles = [pooled_reads_filename]
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)
    graph = finishm_graph.graph
    start_node = finishm_graph.probe_nodes[0]
    start_node_forward = finishm_graph.probe_node_directions[0]
    end_node = finishm_graph.probe_nodes[1]
    end_node_forward = finishm_graph.probe_node_directions[1]

    log.info "Node(s) found that are suitable as initial and terminal nodes in the graph search, respectively: #{start_node.node_id} and #{end_node.node_id}"

    log.info "Removing nodes unconnected to either the start or the end from the graph.."
    original_num_nodes = graph.nodes.length
    original_num_arcs = graph.arcs.length
    filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
    filter.remove_unconnected_nodes(graph, [start_node, end_node])
    log.info "Removed #{original_num_nodes-graph.nodes.length} nodes and #{original_num_arcs-graph.arcs.length} arcs"

    if options[:output_graph_png]
      log.info "Converting assembly to a graphviz PNG"
      viser = Bio::Assembly::ABVisualiser.new
      gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id})
      gv.output :png => options[:output_graph_png], :use => :neato
    end
    if options[:output_graph_svg]
      log.info "Converting assembly to a graphviz SVG"
      viser = Bio::Assembly::ABVisualiser.new
      gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id})
      gv.output :svg => options[:output_graph_svg], :use => :neato
    end
    if options[:output_graph_dot]
      log.info "Converting assembly to a graphviz DOT"
      viser = Bio::Assembly::ABVisualiser.new
      gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id})
      gv.output :dot => options[:output_graph_dot]
    end

    log.info "Searching for trails between the initial and terminal nodes, within the assembly graph"
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    #raise "Untested connection finder below"
    #trails = cartographer.find_all_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
    trails = cartographer.find_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
    log.info "Found #{trails.length} trail(s) between the initial and terminal nodes"

#    log.info "Reading kmer abundances from #{options[:kmer_multiple_abundance_file]}.."
#    kmer_hash = Bio::KmerMultipleAbundanceHash.parse_from_file options[:kmer_multiple_abundance_file]
#    log.info "Finished reading the kmer abundances"

#    if options[:trail_kmer_coverage_file]
#      log.info "Writing out kmer coverages to #{options[:trail_kmer_coverage_file]}.."
#      writer = Bio::AssemblyGraphAlgorithms::KmerCoverageWriter.new
#      io = File.open(options[:trail_kmer_coverage_file],'w')
#      writer.write(io, trails, kmer_hash)
#      log.info "Finished writing"
#    end

#    log.info "Filtering trail(s) based on kmer coverage, requiring each kmer in the path to have a minimum of #{options[:kmer_path_filter_min_coverage]} coverage in patterned reads, except for the #{options[:kmer_path_end_exclusion_length]}bp at the ends"
#    kmer_path_filter = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new
#    thresholds = desired_pattern.collect{|c| c == true ? 1 : 0}
#    log.info "Using thresholds for filtering: #{thresholds}"
#    trails = kmer_path_filter.filter(trails, kmer_hash, thresholds, :exclude_ending_length => options[:kmer_path_end_exclusion_length])
#    log.info "After filtering remained #{trails.length} trails"

    log.debug "Found trails: #{trails.collect{|t| t.to_s}.join("\n")}"

    # TODO: actually output the trails somehow
#    trails.each_with_index do |trail, i|
#      puts ">trail#{i+1}"
#      puts trail.sequence
#    end
  end
end

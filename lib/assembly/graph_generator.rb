require 'bio-velvet'
require 'bio'
require 'pry'

class Bio::FinishM::GraphGenerator
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :velvet_kmer_size => 51,
    :assembly_coverage_cutoff => 3.5,
    }

  def add_options(option_parser, options)
    options.merge!(DEFAULT_OPTIONS)
    option_parser.on("--assembly-kmer NUMBER", "when assembling, use this kmer length [default: #{options[:velvet_kmer_size] }]") do |arg|
      options[:velvet_kmer_size] = arg.to_i
    end
    option_parser.on("--assembly-coverage-cutoff NUMBER", "Require this much coverage in each node, all other nodes are removed [default: #{options[:assembly_coverage_cutoff] }]") do |arg|
      options[:assembly_coverage_cutoff] = arg.to_f
    end
    option_parser.on("--post-assembly-coverage-cutoff NUMBER", "Require this much coverage in each node, implemented after assembly [default: not used]") do |arg|
      options[:post_assembly_coverage_cutoff] = arg.to_f
    end
    option_parser.on("--velvet-directory PATH", "Output assembly intermediate files to this directory [default: use temporary directory, delete afterwards]") do |arg|
      options[:output_assembly_path] = arg
    end
    option_parser.on("--already-assembled-velvet-directory PATH", "If an assembly directory has been specified previously with --velvet-directory, re-use this assembly rather than re-doing the assembly [default: off]") do |arg|
      options[:previous_assembly] = arg
    end
  end

  # Generate a ProbedGraph object, given one or more 'probe sequences'
  # and metagenomic reads. This is a rather large method, but seems to
  # be approximately repeated in different applications of FinishM, so
  # creating it for DRY purposes.
  #
  # probe_sequences: DNA sequences (as String objects whose direction points to the outsides of contigs)
  # read_inputs: a ReadInput object, containing the information to feed to velveth
  #
  # options:
  # :probe_reads: a list of sequence numbers (numbering as per velvet Sequence file)
  # :probe_read_names: a list of sequence names (not IDs) that are probes (convert the names to IDs using the CnyUnifiedSeqNames file)
  # :velvet_kmer_size: kmer
  # :assembly_coverage_cutoff: coverage cutoff for nodes
  # :post_assembly_coverage_cutoff: apply this coverage cutoff to nodes after parsing assembly
  # :output_assembly_path: write assembly to this directory
  # :previous_assembly: a velvet directory from a previous run of the same probe sequences and reads. (Don't re-assemble)
  # :use_textual_sequence_file: by default, a binary sequence file is used. Set this true to get velvet to generate the Sequences file
  # :remove_unconnected_nodes: delete nodes from the graph that are not connected to the probe nodes
  # :graph_search_leash_length: when :remove_unconnected_nodes'ing, use this leash length
  def generate_graph(probe_sequences, read_inputs, options={})
    options[:parse_sequence_file] ||= true
    graph = nil
    read_probing_graph = nil
    finishm_graph = Bio::FinishM::ProbedGraph.new

    velvet_result = nil

    probe_read_ids = nil
    if options[:probe_reads]
      probe_read_ids = options[:probe_reads]
    else
      probe_read_ids = Set.new((1..probe_sequences.length))
    end
    if options[:previous_assembly].nil? #If assembly has not already been carried out
      Tempfile.open('probes.fa') do |tempfile|
        50.times do # Do 50 times to make sure that velvet doesn't throw out parts of the graph that contain this contig
          probe_sequences.each_with_index do |probe, i|
            tempfile.puts ">probe#{i}"
            tempfile.puts probe
          end
        end
        tempfile.close
        singles = read_inputs.fasta_singles
        if singles and !singles.empty?
          read_inputs.fasta_singles = [tempfile.path, singles].flatten
        else
          read_inputs.fasta_singles = [tempfile.path]
        end
        log.debug "Inputting probes into the assembly:\n#{File.open(tempfile.path).read}" if log.debug?

        log.info "Assembling sampled reads with velvet"
        # Bit of a hack, but have to use -short1 as the anchors because then start and end anchors will have node IDs 1,2,... etc.
        use_binary = options[:use_textual_sequence_file] ? '' : '-create_binary'
        velvet_result = Bio::Velvet::Runner.new.velvet(
          options[:velvet_kmer_size],
          "#{read_inputs.velvet_read_arguments} #{use_binary}",
          "-read_trkg yes -cov_cutoff #{options[:assembly_coverage_cutoff] } -tour_bus no",
          :output_assembly_path => options[:output_assembly_path]
          )
        if log.debug?
          log.debug "velveth stdout: #{velvet_result.velveth_stdout}"
          log.debug "velveth stderr: #{velvet_result.velveth_stderr}"
          log.debug "velvetg stdout: #{velvet_result.velvetg_stdout}"
          log.debug "velvetg stderr: #{velvet_result.velvetg_stderr}"
        end
        log.info "Finished running assembly"
        finishm_graph.velvet_result_directory = velvet_result.result_directory
      end
    else
      log.info "Using previous assembly stored in #{options[:previous_assembly] }"
      velvet_result = Bio::Velvet::Result.new
      velvet_result.result_directory = options[:previous_assembly]
      finishm_graph.velvet_result_directory = velvet_result.result_directory
    end

    # Check that the probe reads given are present in the assembly passed here
    sequence_store = parse_velvet_binary_reads(velvet_result.result_directory)
    finishm_graph.velvet_sequences = sequence_store
    if !check_probe_sequences(probe_sequences, sequence_store)
      raise "Probe sequences changed since previous velvet assembly!"
    end

    log.info "Parsing the graph output from velvet"
    opts = {}
    unless options[:parse_all_noded_reads]
      #Ignore parsing reads that are not probes, as we don't care and this just takes up extra computational resources
      opts[:dont_parse_noded_reads] = true
    end
    graph = Bio::Velvet::Graph.parse_from_file(
      File.join(velvet_result.result_directory, 'LastGraph'),
      opts
      )
    log.info "Finished parsing graph: found #{graph.nodes.length} nodes and #{graph.arcs.length} arcs"

    log.info "Beginning parse of graph using velvet's parsing C code.."
    read_probing_graph = Bio::Velvet::Underground::Graph.parse_from_file File.join(velvet_result.result_directory, 'LastGraph')
    log.info "Completed velvet code parsing velvet graph"

    # Find the anchor nodes again
    anchor_sequence_ids = probe_read_ids.to_a.sort
    endings = []
    unless probe_read_ids.empty? and options[:probe_read_names].nil? #don't bother trying to find probes if none exists
      # Convert read names to read IDs if required
      if options[:probe_read_names]
        # Probe reads are given as names, not IDs. What are the corresponding probes then?
        entries = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries_using_grep_hack(
          File.join(velvet_result.result_directory, 'CnyUnifiedSeq.names'),
          options[:probe_read_names]
          )
        if entries.length != options[:probe_read_names].length
          raise "Unexpected number of probes recovered through names: expected #{options[:probe_read_names].length}, found #{entries.length}"
        end
        anchor_sequence_ids = []
        options[:probe_read_names].each do |name| #maintain order of them as they are specified in the original array parameter
          if entries[name].length > 1
            raise "Found >1 sequence named #{name} in the assembly, being conservative and not continuing"
          else
            anchor_sequence_ids.push entries[name][0].read_id
          end
        end
        log.info "Correctly recovered all #{anchor_sequence_ids.length} sequences using their names"
      end

      finder = Bio::FinishM::CProbeNodeFinder.new
      log.info "Finding probe nodes in the assembly"
      c_graph_endings = finder.find_probes(read_probing_graph, anchor_sequence_ids)
      log.debug "Converting probe nodes found in C graph to Ruby analogues and adding to Ruby-parsed graph"
      endings = c_graph_endings.collect do |node_direction_read|
        if node_direction_read.empty?
          # No probe found
          []
        else #found a node.
          #equivalent node
          node = graph.nodes[node_direction_read[0].node_id]
          #equivalent direction
          direction = node_direction_read[1]
          #equivalent noded read
          nr = Bio::Velvet::Graph::NodedRead.new
          #           nr.read_id = read_id
          #           nr.offset_from_start_of_node = row[1].to_i
          #           nr.start_coord = row[2].to_i
          #           nr.direction = current_node_direction
          cnr = node_direction_read[2]
          nr.read_id = cnr.read_id
          nr.offset_from_start_of_node = cnr.offset_from_start_of_node
          nr.start_coord = cnr.start_coord
          nr.direction = direction
          # collect
          [node, direction, nr]
        end
      end
    end
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = endings.collect{|array| array[0]}
    finishm_graph.probe_node_directions = endings.collect{|array| array[1]}
    finishm_graph.probe_node_reads = endings.collect{|array| array[2]}

    # Check to make sure the probe sequences map to nodes in the graph
    if finishm_graph.completely_probed?
      if log.info?
        found_all = true
        num_found = 0
        finishm_graph.probe_nodes.each_with_index do |probe,i|
          if probe.nil?
            found_all = false
            log.debug "Unable to recover probe ##{i+1}, perhaps this will cause problems, but proceding optimistically"
          else
            num_found += 1
          end
        end
        if found_all
          if finishm_graph.probe_nodes.empty?
            log.debug "No probes specified, so didn't find any"
          else
            log.info "Found all anchoring nodes in the graph."
          end
        else
          log.info "Found #{num_found} of #{finishm_graph.probe_nodes.length} anchoring nodes in the graph, ignoring the rest"
        end
      end
    else
      raise "Unable to find all anchor reads from the assembly, cannot continue. This is probably an error with this script, not you. Probes not found: #{finishm_graph.missing_probe_indices.inspect}"
    end

    if options[:post_assembly_coverage_cutoff]
      log.info "Removing nodes with coverage < #{options[:post_assembly_coverage_cutoff] } from graph.."
      original_num_nodes = graph.nodes.length
      original_num_arcs = graph.arcs.length
      filter = Bio::AssemblyGraphAlgorithms::CoverageBasedGraphFilter.new
      filter.remove_low_coverage_nodes(graph,
        options[:post_assembly_coverage_cutoff],
        :whitelisted_sequences => Set.new(anchor_sequence_ids)
      )
      log.info "Removed #{original_num_nodes-graph.nodes.length} nodes and #{original_num_arcs-graph.arcs.length} arcs, leaving #{graph.nodes.length} nodes and #{graph.arcs.length} arcs."
    end

    if options[:remove_unconnected_nodes]
      if options[:graph_search_leash_length]
        log.info "Removing nodes unconnected to probe nodes from the graph using leash #{options[:graph_search_leash_length] }.."
      else
        log.info "Removing nodes unconnected to probe nodes from the graph without using a leash.."
      end
      original_num_nodes = graph.nodes.length
      original_num_arcs = graph.arcs.length
      filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
      filter.remove_unconnected_nodes(
        graph,
        finishm_graph.probe_nodes.reject{|n| n.nil?},
        :leash_length => options[:graph_search_leash_length]
        )
      log.info "Removed #{original_num_nodes-graph.nodes.length} nodes and #{original_num_arcs-graph.arcs.length} arcs, leaving #{graph.nodes.length} nodes and #{graph.arcs.length} arcs."
    end

    return finishm_graph
  end

  # Read in the reads from a velvet result
  def parse_velvet_binary_reads(velvet_result_directory)
    sequences_file_path = File.join velvet_result_directory, 'CnyUnifiedSeq'
    log.info "Reading in the actual sequences of all reads from #{sequences_file_path}"
    sequences = Bio::Velvet::Underground::BinarySequenceStore.new sequences_file_path
    log.info "Read in #{sequences.length} sequences"
    return sequences
  end

  # When re-using an assembly, sometimes need to make
  # sure that the probe sequences used previously are the same
  # as what is given this time. Given am Array of probe sequences
  # and a binary_sequence_file, check the probe sequences are the
  # consistent.
  def check_probe_sequences(probe_sequences, sequence_store)
    return true if probe_sequences.nil?

    probe_sequences.each_with_index do |probe, i|
      log.debug "Checking probe sequence \##{i+1}" if log.debug?
      if sequence_store[i+1].upcase != probe.upcase
        log.error "Probe sequence \##{i+1} has changed - perhaps the wrong velvet assembly directory was specified, or a fresh assembly is required?"
        return false
      end
    end
    log.debug "Presence of #{probe_sequences.length} probe sequences verified"
    return true
  end
end

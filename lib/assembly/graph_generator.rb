require 'bio-velvet'
require 'bio'
require 'pry'
require 'files'

class Bio::FinishM::ProbedGraph
  attr_accessor :probe_nodes, :probe_node_directions, :probe_node_reads, :graph

  attr_accessor :velvet_result_directory

  # Were all the probe recovered through the process?
  def completely_probed?
    !(@probe_nodes.find{|node| node.nil?})
  end

  def missing_probe_indices
    missings = []
    @probe_nodes.each_with_index do |probe, i|
      missings.push(i+1) if probe.nil?
    end
    return missings
  end

  # Make a Bio::Velvet::Graph::OrientedNodeTrail with just one
  # step in it - the node that corresponds to the probe_index
  def initial_path_from_probe(probe_index)
    initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
    node = @probe_nodes[probe_index]
    raise "No node found for probe #{probe_index}" if node.nil?
    direction = @probe_node_directions[probe_index]

    way = direction ?
    Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST :
      Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
    initial_path.add_node node, way
    return initial_path
  end

  # Return a Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode
  # corresponding to the index of the probe and its direction
  def velvet_oriented_node(probe_index)
    node = @probe_nodes[probe_index]
    if node.nil?
      return nil
    else
      return initial_path_from_probe(probe_index)[0]
    end
  end

  # The leash is the number of base pairs from the start of the probe,
  # but the path finding algorithm simply uses the combined length of all
  # the nodes without reference to the actual probe sequence. So if the
  # probe is near the end of a long node, then path finding may fail.
  # So adjust the leash length to account for this (or keep the nil
  # if the starting_leash_length is nil)
  def adjusted_leash_length(probe_index, starting_leash_length)
    return nil if starting_leash_length.nil?

    read = @probe_node_reads[probe_index]
    return read.offset_from_start_of_node+starting_leash_length
  end
end

# A class representing reads or sets of reads to be assembled
class Bio::FinishM::ReadInput
  READ_INPUT_SYMBOLS = [
    :fasta_singles, :fastq_singles, :fasta_singles_gz, :fastq_singles_gz,
    #:interleaved_fasta, :interleaved_fastq, :interleaved_fasta_gz, :interleaved_fastq_gz,
    ]
  READ_INPUT_SYMBOLS.each do |sym|
    attr_accessor sym
  end

  # Given an OptionParser, add options to it, which parse out read-related options
  def add_options(option_parser, options)
    {
      '--fasta' => :fasta_singles,
      '--fastq' => :fastq_singles,
      '--fasta-gz' => :fasta_singles_gz,
      '--fastq-gz' => :fastq_singles_gz,
      '--interleaved-fasta' => :interleaved_fasta, #paired-end requires more fiddling with the velvet source to remove tip clipping
      '--interleaved-fastq' => :interleaved_fastq,
      '--interleaved-fasta-gz' => :interleaved_fasta_gz,
      '--interleaved-fastq-gz' => :interleaved_fastq_gz,
      }.each do |flag, sym|
        option_parser.on("#{flag} PATH", Array, "One or more paths to reads, comma separated") do |arg|
          options[sym] = arg
        end
      end
  end

  # Require at least 1 set of reads to be given, of any type
  def validate_options(options, argv)
    return nil if options[:previous_assembly] #bit of a hack, but hey
    READ_INPUT_SYMBOLS.each do |sym|
      return nil if options[sym]
    end
    return "No definition of reads for assembly was found"
  end

  # Parse options from options hash into instance variables for this object
  def parse_options(options)
    READ_INPUT_SYMBOLS.each do |sym|
      send("#{sym}=",options[sym]) if options[sym]
    end
  end

  # Like velvet_read_arguments but arguments are
  # returned in an Array, not a String
  def velvet_read_arguments_array
    args = []
    #Put paired sequences first in the hash (in Ruby, hashes are ordered) so that if they are paired, then odd numbered sequences
    # are read1 and even numbered sequences are read2
    {
#       :interleaved_fasta => %w(-fasta -shortPaired),
#       :interleaved_fastq => %w(-fastq -shortPaired),
#       :interleaved_fasta_gz => %w(-fasta.gz -shortPaired),
#       :interleaved_fastq_gz => %w(-fastq.gz -shortPaired),
      :fasta_singles => %w(-fasta -short),
      :fastq_singles => %w(-fastq -short),
      :fasta_singles_gz => %w(-fasta.gz -short),
      :fastq_singles_gz => %w(-fastq.gz -short),
      }.each do |sym, velvet_flags|
        paths = send(sym)
        unless paths.nil? or paths.empty?
          args.push velvet_flags
          paths.each do |path|
            args.push path
          end
        end
      end
    return args.flatten
  end

  # Output a string to be used on the command line with velvet
  def velvet_read_arguments
    velvet_read_arguments_array.collect{|s| s.inspect}.join(' ')
  end
end

class Bio::FinishM::GraphGenerator
  include Bio::FinishM::Logging

  def add_options(option_parser, options)
    options.merge!({
      :velvet_kmer_size => 51,
      :assembly_coverage_cutoff => 3.5,
      })
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
    option_parser.on("--already-assembled-velvet-directory PATH", "Skip until after assembly in this process, and start from this assembly directory created during a previous run of this script [default: off]") do |arg|
      options[:previous_assembly] = arg
    end
    #         option_parser.on("--serialize-velvet-graph FILE", "So that the velvet graph does not have to be reparsed, serialise the parsed object for later use in this file [default: off]") do |arg|
    #           options[:serialize_parsed_graph_file] = arg
    #         end
    #         option_parser.on("--already-serialized-velvet-graph FILE", "Restore the parsed velvet graph from this file [default: off]") do |arg|
    #           options[:previously_serialized_parsed_graph_file] = arg
    #         end
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
  # :velvet_kmer_size: kmer
  # :assembly_coverage_cutoff: coverage cutoff for nodes
  # :post_assembly_coverage_cutoff: apply this coverage cutoff to nodes after parsing assembly
  # :output_assembly_path: write assembly to this directory
  # :previous_assembly: a velvet directory from a previous run of the same probe sequences and reads. (Don't re-assemble)
  # :use_textual_sequence_file: by default, a binary sequence file is used. Set this true to get velvet to generate the Sequences file
  # :remove_unconnected_nodes: delete nodes from the graph that are not connected to the probe nodes
  # :graph_search_leash_length: when :remove_unconnected_nodes'ing, use this leash length
  # :serialize_parsed_graph_file: after assembly, parse the graph in, and serialize the ruby object for later. Value of this option is the path to the save file.
  # :previously_serialized_parsed_graph_file: read in a previously serialized graph file, and continue from there
  # :parse_all_noded_reads: parse NodedRead objects for all nodes (default: only worry about the nodes containing the probe reads)
  def generate_graph(probe_sequences, read_inputs, options={})
    graph = nil
    read_probing_graph = nil
    finishm_graph = Bio::FinishM::ProbedGraph.new

    if options[:previously_serialized_parsed_graph_file].nil?
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

          # Velvet directory is either
          velvet_directory_for_underground = options[:output_assembly_path]
          velvet_directory_for_underground ||= Files.create.root

          Bio::Velvet::Underground::Runner.run(
            options[:velvet_kmer_size],
            read_inputs.velvet_read_arguments_array + [use_binary],
            %w(-read_trkg yes -cov_cutoff)+[options[:assembly_coverage_cutoff]] + %w(-tour_bus no -clip_tips no),
            :velvet_directory => velvet_directory_for_underground#options[:output_assembly_path]
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
        log.info "Using previous assembly stored at #{options[:previous_assembly] }"
        velvet_result = Bio::Velvet::Result.new
        velvet_result.result_directory = options[:previous_assembly]
        finishm_graph.velvet_result_directory = velvet_result.result_directory
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

      if options[:serialize_parsed_graph_file]
        log.info "Storing a binary version of the graph file for later use at #{options[:serialize_parsed_graph_file] }"
        File.open(options[:serialize_parsed_graph_file],'wb') do |f|
          f.print Marshal.dump(graph)
        end
        log.info "Stored a binary representation of the velvet graph at #{options[:serialize_parsed_graph_file] }"
      end
    else
      log.info "Restoring graph file from #{options[:previously_serialized_parsed_graph_file] }.."
      graph = Marshal.load(File.open(options[:previously_serialized_parsed_graph_file]))
      log.info "Restoration complete"
    end


    # Find the anchor nodes again
    anchor_sequence_ids = probe_read_ids.to_a.sort
    endings = []
    unless probe_read_ids.empty? #don't bother trying to find probes if none exists
      finder = Bio::AssemblyGraphAlgorithms::NodeFinder.new
      log.info "Finding probe nodes in the assembly"
      endings = finder.find_unique_nodes_with_sequence_ids(read_probing_graph, anchor_sequence_ids)
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
end

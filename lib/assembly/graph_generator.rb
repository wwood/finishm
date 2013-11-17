require 'bio-velvet'
require 'bio'

module Bio
  module FinishM
    class ProbedGraph
      attr_accessor :probe_nodes, :probe_node_directions, :graph

      # Were all the probe recovered through the process?
      def completely_probed?
        !(@probe_nodes.find{|node| node.nil?}
      end

      def missing_probe_indices
        missings = []
        @probe_nodes.each_with_index do |probe, i|
          missings.push(i+1) if probe.nil?
        end
        return missings
      end
    end

    # A class representing reads or sets of reads to be assembled
    class ReadInput
      attr_accessor :fasta_singles, :fastq_singles, :fasta_singles_gz, :fastq_singles_gz

      # Given an OptionParser, add options to it, which parse out read-related options
      def add_options(option_parser, options)
        {
          '--fasta' => :fasta_singles,
          '--fastq' => :fastq_singles,
          '--fasta-gz' => :fasta_singles_gz,
          '--fastq-gz' => :fastq_singles_gz,
        }.each do |flag, sym|
          option_parser.on("#{flag} PATH") do |arg|
            send("#{sym.to_s}||=".to_sym, []) #initialise the instance variable as an empty array

            # add entries to the instance variable array
            attr_symbol = "#{sym.to_s}||=".to_sym
            arg.split(/[^\\] /).each do |filename|
              filename.strip
              send(attr_symbol).push filename
            end
          end
        end
      end

      def velvet_read_arguments
        readset_index = 1
        args = ''
        {
          :fasta_singles => '-fasta',
          :fastq_singles => '-fastq',
          :fasta_singles_gz => '-fasta.gz',
          :fastq_singles_gz => '-fastq.gz'
        }.each do |sym, velvet_flag|
          paths = send(sym)
          unless paths.nil? or paths.empty?
            args += " #{velvet_flag}"
            paths.each do |path|
              args += "-short#{readset_index} #{path}"
              readset_index += 1
            end
          end
        end
        return args
      end
    end

    class GraphGenerator
      include Bio::FinishM::Logging

      def add_options(option_parser, options)
        option_parser.on("--assembly-kmer NUMBER", "when assembling, use this kmer length [default: #{options[:velvet_kmer_size]}]") do |arg|
          options[:velvet_kmer_size] = arg.to_i
        end
        option_parser.on("--assembly-coverage-cutoff NUMBER", "Require this much coverage in each node, all other nodes are removed [default: #{options[:assembly_coverage_cutoff]}]") do |arg|
          options[:assembly_coverage_cutoff] = arg.to_f
        end
        option_parser.on("--velvet-directory PATH", "Output assembly intermediate files to this directory [default: off]") do |arg|
          options[:output_assembly_path] = arg
        end
        option_parser.on("--already-assembled-velvet-directory PATH", "Skip until after assembly in this process, and start from this assembly directory created during a previous run of this script [default: off]") do |arg|
          options[:previous_assembly] = arg
        end
        option_parser.on("--serialize-velvet-graph FILE", "So that the velvet graph does not have to be reparsed, serialise the parsed object for later use in this file [default: off]") do |arg|
          options[:serialize_parsed_graph_file] = arg
        end
        option_parser.on("--already-serialized-velvet-graph FILE", "Restore the parsed velvet graph from this file [default: off]") do |arg|
          options[:previously_serialized_parsed_graph_file] = arg
        end
      end

      # Generate a ProbedGraph object, given one or more 'probe sequences'
      # and metagenomic reads. This is a rather large method, but seems to
      # be approximately repeated in different applications of FinishM, sor
      # creating it for DRY purposes.
      #
      # probe_sequences: DNA sequences (as String objects whose direction points to the outsides of contigs)
      # read_inputs: a ReadInput object, containing the information to feed to velveth
      #
      # options:
      # :velvet_kmer_size: kmer
      # :assembly_coverage_cutoff: coverage cutoff for nodes
      # :output_assembly_path: write assembly to this directory
      # :previous_assembly: a velvet directory from a previous run of the same probe sequences and reads. (Don't re-assemble)
      # :serialize_parsed_graph_file: after assembly, parse the graph in, and serialize the ruby object for later. Value of this option is the path to the save file.
      # :previously_serialized_parsed_graph_file: read in a previously serialized graph file, and continue from there
      def generate_graph(probe_sequences, read_inputs, options)
        graph = nil
        if options[:previously_serialized_parsed_graph_file].nil?
          velvet_result = nil
          if options[:previous_assembly].nil? #If assembly has not already been carried out
            Tempfile.open('probes.fa') do |tempfile|
              probe_sequences.each_with_index do |probe, i|
                tempfile.puts ">probe#{i}"
                tempfile.puts probe
              end
              tempfile.close
              log.debug "Inputting probes into the assembly:\n#{File.open(tempfile.path).read}" if log.debug?

              log.info "Assembling sampled reads with velvet"
              # Bit of a hack, but have to use -short1 as the anchors because then start and end anchors will have node IDs 1,2,... etc.
              velvet_result = Bio::Velvet::Runner.new.velvet(
                options[:velvet_kmer_size],
                read_inputs.velvet_read_arguments,
                "-read_trkg yes -cov_cutoff #{options[:assembly_coverage_cutoff]}",
                :output_assembly_path => options[:output_assembly_path]
              )
              if log.debug?
                log.debug "velveth stdout: #{velvet_result.velveth_stdout}"
                log.debug "velveth stderr: #{velvet_result.velveth_stderr}"
                log.debug "velvetg stdout: #{velvet_result.velvetg_stdout}"
                log.debug "velvetg stderr: #{velvet_result.velvetg_stderr}"
              end
              log.info "Finished running assembly"
            end
          else
            log.info "Using previous assembly stored at #{options[:previous_assembly]}"
            velvet_result = Bio::Velvet::Result.new
            velvet_result.result_directory = options[:previous_assembly]
          end

          log.info "Parsing the graph output from velvet"
          graph = Bio::Velvet::Graph.parse_from_file(File.join velvet_result.result_directory, 'LastGraph')
          log.info "Finished parsing graph: found #{graph.nodes.length} nodes and #{graph.arcs.length} arcs"

          if options[:serialize_parsed_graph_file]
            log.info "Storing a binary version of the graph file for later use at #{options[:serialize_parsed_graph_file]}"
            File.open(options[:serialize_parsed_graph_file],'wb') do |f|
              f.print Marshal.dump(graph)
            end
            log.info "Stored a binary representation of the velvet graph at #{options[:serialize_parsed_graph_file]}"
          end
        else
          log.info "Restoring graph file from #{options[:previously_serialized_parsed_graph_file]}.."
          graph = Marshal.load(File.open(options[:previously_serialized_parsed_graph_file]))
          log.info "Restoration complete"
        end


        # Find the anchor nodes again
        finder = Bio::AssemblyGraphAlgorithms::NodeFinder.new
        log.info "Finding nodes representing the end of the each contig"
        anchor_sequence_ids = (1..probe_sequences.length)
        endings = finder.find_unique_nodes_with_sequence_ids(graph, anchor_sequence_ids)
        finishm_graph = Bio::FinishM::ProbedGraph.new
        endings.each do |array|
          node = array[0]
          direction = array[1]

          finishm_graph.probe_nodes ||= []
          finishm_graph.probe_nodes.push node
          finishm_graph.probe_node_directions ||= []
          finishm_graph.probe_node_directions.push direction
        end

        # Check to make sure the probe sequences map to nodes in the graph
        if finishm_graph.completely_probed?
          probe_descriptions = finishm_graph.probe_nodes.collect do |probe,i|
            "#{probe.node_id}/#{finishm_graph.probe_node_directions[i]}"
          end.join(', ')
          log.info "Found all anchoring nodes in the graph: #{probe_descriptions}"
        else
          raise "Unable to find both anchor reads from the assembly, cannot continue. This is probably an error with this script, not you. Probes not found: #{finishm_graph.missing_probe_indices.inspect}"
        end

        log.info "Removing nodes unconnected to either the start or the end from the graph.."
        original_num_nodes = graph.nodes.length
        original_num_arcs = graph.arcs.length
        filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
        filter.remove_unconnected_nodes(graph, finishm_graph.probe_nodes)
        log.info "Removed #{original_num_nodes-graph.nodes.length} nodes and #{original_num_arcs-graph.arcs.length} arcs, leaving #{graph.nodes.length} nodes and #{graph.arcs.length} arcs."

        return finishm_graph
      end
    end
  end
end

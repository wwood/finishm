#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-velvet'
require 'tempfile'
require 'pp'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'finishm'
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'priner'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :velvet_kmer_size => 73,#TODO: these options should be exposed to the user, and perhaps not guessed at
  :velvetg_arguments => '-read_trkg yes',# -exp_cov 41 -cov_cutoff 12.0973243610491', #hack
  :contig_end_length => 300,
  :output_assembly_path => 'velvetAssembly',
  :graph_search_leash_length => 3000,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} --reads <read_file> --contigs <contigs_file>

    Takes a set of reads and a set of contigs. Then it runs an assembly based on those reads,
    and tries to fill in possible gaps between the contigs. There may be multiple ways
    to join two contig ends together - in this that multiple cases are reported. \n\n"


  opts.on("--reads FILE", "gzipped fastq file of reads to perform the re-assembly with [required]") do |arg|
    options[:reads_file] = arg
  end
  opts.on("--contigs FILE", "fasta file of contigs to be joined together [required]") do |arg|
    options[:contigs_file] = arg
  end

  opts.separator "\nOptional arguments:\n\n"
  opts.on("--output-trails-fasta PATH", "Output found paths to this file in fasta format [default: off]") do |arg|
    options[:overall_trail_output_fasta_file] = arg
  end
  opts.on("--already-assembled-velvet-directory PATH", "Skip until after assembly in this process, and start from this assembly directory created during a previous run of this script [default: off]") do |arg|
    options[:previous_assembly] = arg
  end
  opts.on("--serialize-velvet-graph FILE", "So that the velvet graph does not have to be reparsed, serialise the parsed object for later use in this file [default: off]") do |arg|
    options[:serialize_parsed_graph_file] = arg
  end
  opts.on("--already-serialized-velvet-graph FILE", "Restore the parsed velvet graph from this file [default: off]") do |arg|
    options[:previously_serialized_parsed_graph_file] = arg
  end


  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:reads_file].nil? or options[:contigs_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
Bio::Log::LoggerPlus.new 'bio-velvet'; Bio::Log::CLI.configure 'bio-velvet'

# Extract contig ends from each of the input contigs, so that the contig ends can be found in the
# assembly graph structure.
contig_ends = []
velvet_sequence_id_to_contig_end = {}
contig_lengths = {}
class ContigEnd
  attr_accessor :sequence, :start_or_end, :contig_name, :velvet_sequence_id
end
velvet_read_index = 1
Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
  contig_lengths[seq.definition] = seq.seq.length
  if seq.seq.length < options[:contig_end_length]
    log.warn "Contig #{seq.definition} is shorter than the end length used to anchor the contig in the assembly. This is not ideal but may be ok."
    #TODO: fix this - should be counting from the middle. Should I just ignore those ones?
  end
  # Add the start of the contig
  contig_end = ContigEnd.new
  contig_end.start_or_end = :start
  contig_end.sequence = Bio::Sequence::NA.new(seq.seq[0...options[:contig_end_length]]).reverse_complement.to_s
  contig_end.contig_name = seq.definition
  velvet_sequence_id_to_contig_end[velvet_read_index] = contig_end
  contig_end.velvet_sequence_id = velvet_read_index; velvet_read_index += 1
  contig_ends.push contig_end


  # Add the back of the contig
  contig_end = ContigEnd.new
  contig_end.start_or_end = :end
  s = seq.seq
  contig_end.sequence = s[s.length-options[:contig_end_length]...s.length]
  contig_end.contig_name = seq.definition
  velvet_sequence_id_to_contig_end[velvet_read_index] = contig_end
  contig_end.velvet_sequence_id = velvet_read_index; velvet_read_index += 1
  contig_ends.push contig_end
end
log.info "Parsed in #{contig_ends.length} contig ends from the two sides of each input contig"


graph = nil
if options[:previously_serialized_parsed_graph_file].nil?
  velvet_result = nil
  if options[:previous_assembly].nil? #If assembly has not already been carried out
    Tempfile.open('anchors.fa') do |tempfile|
      contig_ends.each do |contig_end|
        tempfile.puts ">anchor#{contig_end.velvet_sequence_id}"
        tempfile.puts contig_end.sequence
      end

      log.info "Assembling sampled reads with velvet"
      # Bit of a hack, but have to use -short1 as the anchors because then start and end anchors will have node IDs 1,2,... etc.
      velvet_result = Bio::Velvet::Runner.new.velvet(
        options[:velvet_kmer_size],
        "-short #{tempfile.path} -short2 -fastq.gz #{options[:reads_file]}",
        options[:velvetg_arguments],
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
  graph = Bio::Velvet::Graph.parse_from_file(File.join velvet_result.result_directory, 'Graph2')
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



# Find the anchoring nodes for each of the contig ends
finder = Bio::AssemblyGraphAlgorithms::NodeFinder.new
log.info "Finding node representing the end of the each contig"
i = 1
anchor_sequence_ids = contig_ends.collect{|c| c.velvet_sequence_id}
anchoring_nodes_and_directions = finder.find_unique_nodes_with_sequence_ids(graph, anchor_sequence_ids)
num_anchors_found = anchoring_nodes_and_directions.reject{|s,e| e[0].nil?}.length
anchoring_node_id_to_contig_end = {}
anchoring_nodes_and_directions.each do |seq_id, node_and_direction|
  next if node_and_direction[0].nil? #skip when there is no node found in the graph for this contig end
  anchoring_node_id_to_contig_end[node_and_direction[0].node_id] = velvet_sequence_id_to_contig_end[seq_id]
end
log.info "Found anchoring nodes for #{num_anchors_found} out of #{contig_ends.length} contig ends"

log.info "Searching for trails between the nodes within the assembly graph"
cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
trail_sets = cartographer.find_trails_between_node_set(graph, anchoring_nodes_and_directions.values.reject{|v| v[0].nil?}, options[:graph_search_leash_length])
log.info "Found #{trail_sets.reduce(0){|s,set|s+=set.length}} trail(s) in total"

node_id_to_contig_description = {}
anchoring_nodes_and_directions.each do |seq_id, pair|
  next if pair.empty? #When no nodes were found
  node_id = pair[0].node_id
  node_id_to_contig_description[node_id] = velvet_sequence_id_to_contig_end[seq_id]
end
contig_end_id_to_partners = {}
# Tabulate all the partners each way (complete the previously triangular matrix)
trail_sets.each do |trail_set|
  trail_set.each do |trail|
    start_id = trail.first.node.node_id
    end_id = trail.last.node.node_id
    contig_end_id_to_partners[start_id] ||= []
    contig_end_id_to_partners[start_id].push node_id_to_contig_description[end_id]
    contig_end_id_to_partners[end_id] ||= []
    contig_end_id_to_partners[end_id].push node_id_to_contig_description[start_id]
  end
end

puts %w(contig_end_id contig_name contig_length connections).join "\t"
trail_sets.each_with_index do |trail_set, i|
  partner_contig_ends = contig_end_id_to_partners[contig_ends[i].velvet_sequence_id]
  partner_contig_ends ||= []
  # Each contig has 2 trail sets associated with it - one for the start and one for the end
  puts [
    contig_ends[i].velvet_sequence_id,
    contig_ends[i].contig_name,
    contig_lengths[contig_ends[i].contig_name],
    partner_contig_ends.collect{|c| c.velvet_sequence_id}.sort.join(',')
  ].join("\t")
end

if options[:overall_trail_output_fasta_file]
  File.open(options[:overall_trail_output_fasta_file],'w') do |outfile|
    trail_sets.each do |trail_set|
      trail_set.each do |trail|
        begin
          trail_sequence = trail.sequence #Get the trail sequence first as this may not be possible.

          start_id = trail.first.node.node_id
          end_id = trail.last.node.node_id
          start_contig_end = anchoring_node_id_to_contig_end[start_id]
          end_contig_end = anchoring_node_id_to_contig_end[end_id]
          outfile.print '>'
          outfile.print start_contig_end.contig_name
          outfile.print '_'
          outfile.print start_contig_end.start_or_end
          outfile.print ':'
          outfile.print end_contig_end.contig_name
          outfile.print '_'
          outfile.puts end_contig_end.start_or_end

          outfile.puts trail_sequence
        rescue Bio::Velvet::NotImplementedException => e
          log.warn "Problem getting sequence of found trail #{trail.to_s}, skipping this trail: #{e.to_s}"
        end
      end
    end
  end
end



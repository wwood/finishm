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
  :velvet_kmer_size => 43,#TODO: these options should be exposed to the user, and perhaps not guessed at
  :contig_end_length => 200,
  :output_assembly_path => nil,#'velvetAssembly',
  :graph_search_leash_length => 3000,
  :assembly_coverage_cutoff => 1.5,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} --reads <read_file> --contig <contig_file>

    Takes a set of reads and a contig that contains gap characters. Then it tries to fill in
    these N characters. It is possible that there is multiple ways to close the gap - in that case
    each is reported. \n\n"


  opts.on("--reads FILE", "gzipped fastq file of reads to perform the gap closing with [required]") do |arg|
    options[:reads_file] = arg
  end
  opts.on("--contig FILE", "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
    options[:contig_file] = arg
  end
  opts.on("--output-trails-fasta PATH", "Output found paths to this file in fasta format [default: off]") do |arg|
    options[:overall_trail_output_fasta_file] = arg
  end

  opts.separator "\nOptional arguments:\n\n"
  opts.on("--overhang NUM", "Start assembling this far from the gap [default: #{options[:contig_end_length]}]") do |arg|
    options[:contig_end_length] = arg.to_i
  end
  opts.on("--start OFFSET", "Start trying to fill from this position in the contig, requires --stop [default: found from position of Ns}]") do |arg|
    options[:start_offset] = arg.to_i-1
  end
  opts.on("--stop OFFSET", "Start trying to fill to this position in the contig, requires --start [default: found from position of Ns}]") do |arg|
    options[:end_offset] = arg.to_i-1
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
  opts.on("--velvet-kmer KMER", "kmer size to use with velvet [default: #{options[:velvet_kmer_size]}]") do |arg|
    options[:velvet_kmer_size] = arg.to_i
  end

  opts.separator "\nDebug-related options:\n\n"



  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:reads_file].nil? or options[:contig_file].nil? or options[:overall_trail_output_fasta_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
Bio::Log::LoggerPlus.new 'bio-velvet'; Bio::Log::CLI.configure 'bio-velvet'

# Find where the Ns are
n_region_start = nil
n_region_end = nil
sequence = nil
Bio::FlatFile.foreach(options[:contig_file]) do |seq|
  if sequence
    raise Exception, "Sorry, this script can only handle single sequences to be gap filled at the moment"
  end

  sequence = seq.seq

  if options[:start_offset] and options[:end_offset]
    log.info "Trying to gap fill from #{options[:start_offset]+1} to #{options[:end_offset]+1}"
    n_region_start = options[:start_offset]
    n_region_end = options[:end_offset]
  else
    log.info "Determining where to fill from the presence of Ns"

    matches = sequence.match(/(N+)/i)
    if !matches
      raise "Unable to find any gaps in the input sequence. That was a bit too easy.."
    end
    n_region_start = matches.offset(0)[0]
    n_region_end = n_region_start + matches[1].length
    log.info "Detected a gap between #{n_region_start} and #{n_region_end}"
  end

  # Check to make sure we are sufficiently distant from the ends
  if n_region_start < options[:contig_end_length] or
    sequence.length - n_region_end < options[:contig_end_length]
    raise "The gap is too close to the end of the contig, sorry"
  end
end

fwd2 = Bio::Sequence::NA.new(sequence[n_region_end..(n_region_end+options[:contig_end_length])])
probe_sequences = [
  sequence[n_region_start-options[:contig_end_length]-1...n_region_start],
  fwd2.reverse_complement.to_s
]


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
  gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id, :digraph => false})
  gv.output :dot => options[:output_graph_dot]
end



log.info "Searching for trails between the nodes within the assembly graph"
cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
trails = cartographer.find_all_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
log.info "Found #{trails.length} trail(s) in total"


log.debug "Outputing trail sequences"
File.open(options[:overall_trail_output_fasta_file],'w') do |f|
  trails.each_with_index do |trail, i|
    f.puts ">trail#{i+1}"
    f.puts trail.sequence
  end
end


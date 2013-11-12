#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'tempfile'
require 'pp'
require 'systemu'
require 'bio-velvet'
require 'set'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'finishm'
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'priner'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :min_leftover_length => false,
  :kmer_coverage_target => 1,
  :velvet_kmer_size => 155,
  :contig_end_length => 300,
  :graph_search_leash_length => 20000,
  :reads_to_assemble => nil,
  :assembly_coverage_cutoff => 1.5,
  :kmer_path_filter_min_coverage => 1,
  :kmer_path_end_exclusion_length => 50,
  :trail_kmer_coverage_file => 'trail_coverages.csv'
}

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

o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} --trails <trail(s).fa> --kmer-abundances <abundances.csv>

    Given an input kmer then abundances space separated file, and a threshold, print out how many kmers are unique to different subsets of columns\n\n"

  opts.on("--trails FASTA_FILE", "fasta file of trails to be tested [required]") do |arg|
    options[:trails_fasta_file] = arg
  end
  opts.on("--kmer-abundances FILE", "kmer multiple abundance file [required]") do |arg|
    options[:kmer_multiple_abundance_file] = arg
  end

  opts.separator "\nOptional arguments:\n\n"

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:upper_threshold].nil? or options[:lower_threshold].nil? or options[:pattern].nil? or options[:kmer_multiple_abundance_file].nil? or options[:reads_files].nil?
  pp options
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

trails_fasta_file = ARGV[0]
# read in fasta file of trails
fasta_seqs = {}
Bio::FlatFile.foreach(options[:trails_fasta_file]) do |e|
  name = e.definition

  e.definition

log.info "Searching for trails between the initial and terminal nodes, within the assembly graph"
cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#raise "Untested connection finder below"
#trails = cartographer.find_all_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
trails = cartographer.find_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
log.info "Found #{trails.length} trail(s) between the initial and terminal nodes"

log.info "Reading kmer abundances from #{options[:kmer_multiple_abundance_file]}.."
kmer_hash = Bio::KmerMultipleAbundanceHash.parse_from_file options[:kmer_multiple_abundance_file]
log.info "Finished reading the kmer abundances"

if options[:trail_kmer_coverage_file]
  log.info "Writing out kmer coverages to #{options[:trail_kmer_coverage_file]}.."
  writer = Bio::AssemblyGraphAlgorithms::KmerCoverageWriter.new
  io = File.open(options[:trail_kmer_coverage_file],'w')
  writer.write(io, trails, kmer_hash)
  log.info "Finished writing"
end

log.info "Filtering trail(s) based on kmer coverage, requiring each kmer in the path to have a minimum of #{options[:kmer_path_filter_min_coverage]} coverage in patterned reads, except for the #{options[:kmer_path_end_exclusion_length]}bp at the ends"
kmer_path_filter = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new
thresholds = desired_pattern.collect{|c| c == true ? 1 : 0}
log.info "Using thresholds for filtering: #{thresholds}"
trails = kmer_path_filter.filter(trails, kmer_hash, thresholds, :exclude_ending_length => options[:kmer_path_end_exclusion_length])
log.info "After filtering remained #{trails.length} trails"

log.debug "Found trails: #{trails.collect{|t| t.to_s}.join("\n")}"

trails.each_with_index do |trail, i|
  puts ">trail#{i+1}"
  puts trail.sequence
end

#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'pp'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'finishm'
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'priner'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}

o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} --trails <trail(s) fasta> --kmer-abundances <abundances.csv>

    Given an input kmer set of sequences then abundances space separated file, and a threshold, print out how many kmers are unique to different subsets of columns\n\n"

  opts.on("--trails FASTA_FILE", "fasta file of trail(s) to be tested [required]") do |arg|
    options[:trails_fasta_file] = arg
  end
  opts.on("--kmer-abundances FILE", "kmer multiple abundance file [required]") do |arg|
    options[:kmer_multiple_abundance_file] = arg
  end

  #opts.separator "\nOptional arguments:\n\n"
  opts.on("--output-kmer-coverages FILE", "output kmer coverages across each library across the contigs default: don't output]") do |arg|
    options[:output_trail_kmer_coverage_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:trails_fasta_file].nil? or options[:kmer_multiple_abundance_file].nil? or options[:output_trail_kmer_coverage_file].nil?
  $stderr.puts "Options not correctly specified. Found:"
  pp options
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

# read in fasta file of trails
log.info "Reading trail sequences from #{options[:trails_fasta_file]}"
fasta_seqs = {}
Bio::FlatFile.foreach(options[:trails_fasta_file]) do |e|
  name = e.definition
  seq = e.seq.seq
  fasta_seqs[name]=seq
end

log.info "Reading kmer abundances from #{options[:kmer_multiple_abundance_file]}.."
kmer_hash = Bio::KmerMultipleAbundanceHash.parse_from_file options[:kmer_multiple_abundance_file]
log.info "Finished reading the kmer abundances"

if options[:output_trail_kmer_coverage_file]
  log.info "Writing out kmer coverages to #{options[:output_trail_kmer_coverage_file]}.."
  writer = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new
  io = File.open(options[:output_trail_kmer_coverage_file],'w')
  fasta_seqs.each do |name, seq|
    log.debug "Writing coverages for #{name}"
    writer.write_depths(io, name, seq, kmer_hash)
  end
  log.info "Finished writing"
end

#log.info "Filtering trail(s) based on kmer coverage, requiring each kmer in the path to have a minimum of #{options[:kmer_path_filter_min_coverage]} coverage in patterned reads, except for the #{options[:kmer_path_end_exclusion_length]}bp at the ends"
#kmer_path_filter = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new
#thresholds = desired_pattern.collect{|c| c == true ? 1 : 0}
#log.info "Using thresholds for filtering: #{thresholds}"
#trails = kmer_path_filter.filter(trails, kmer_hash, thresholds, :exclude_ending_length => options[:kmer_path_end_exclusion_length])
#log.info "After filtering remained #{trails.length} trails"

#trails.each_with_index do |trail, i|
#  puts ">trail#{i+1}"
#  puts trail.sequence
#end

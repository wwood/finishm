#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'bio'
require 'tempfile'
require 'systemu'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'kmer_abundance_pattern'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :number_of_kmers => 100,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} -f <scaffolds_fasta> -k <kmer_abundance_file>

    Take a fasta file of contigs, and a multiple kmer count file. Output the patterns each contig end shows up.n\n"


  opts.on("-f FASTA_FILE", "Fasta file containing multiple sequences that we are attempting to scaffold together [required]") do |arg|
    options[:fasta_file] = arg
  end
  opts.on("-k KMER_FILE", "kmer frequencies [required]") do |arg|
    options[:kmer_file] = arg
  end
  opts.on("--kmer KMER_SIZE", "kmer length [required]") do |arg|
    options[:kmer_size] = arg.to_i
  end
  opts.on("--upper-threshold ARG", "kmer frequency cutoff to saying 'present' [required]") do |arg|
    options[:upper_threshold] = arg.to_i
  end
  opts.on("--lower-threshold ARG", "kmer frequency cutoff to saying 'not present' [required]") do |arg|
    options[:lower_threshold] = arg.to_i
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:fasta_file].nil? or options[:kmer_file].nil? or options[:upper_threshold].nil? or options[:lower_threshold].nil? or options[:kmer_size].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

get_kmer_abundances_from_kmers = lambda do |kmers|
  # fgrep the kmer abundance file for these particular ones
  patterns = {}
  Tempfile.open('kemrs') do |tempfile|
    tempfile.puts kmers.join "\n"
    tempfile.close

    # for each of the kmers that come back, output their pattern in the kmer abundance file
    grep_cmd = "fgrep -f #{tempfile.path} #{options[:kmer_file].inspect}"
    log.debug "Running cmd with #{kmers.length} kmers: #{grep_cmd}"
    status, stdout, stderr = systemu grep_cmd
    raise stderr if stderr != ''
    raise unless status.exitstatus == 0
    num_kmers = stdout.split("\n").length
    log.debug "Finished grepping for kmers, found #{num_kmers} kmers"



    stdout.each_line do |line|
      CSV.parse(line, :col_sep => ' ') do |row|
        pattern = KmerAbundancePattern.new
        pattern.parse_from_kmer_abundance row[1...row.length].collect{|a| a.to_f}, options[:lower_threshold], options[:upper_threshold]
        rep = pattern.binary_string
        patterns[rep] ||= 0
        patterns[rep] += 1
      end
    end
  end
  patterns.sort{|a,b| b[1]<=>a[1]}.collect{|a| a.join(',')}.join("\t")
end

# For each of the sequences in the fasta file
Bio::FlatFile.foreach(Bio::FastaFormat, options[:fasta_file]) do |seq|
  # Extract the first 100 kmers
  kmers = []
  i = 0
  seq.to_biosequence.window_search(options[:kmer_size]) do |s|
    kmers.push s.seq
    i += 1
    break if i>options[:number_of_kmers]
  end
  # output to a tempfile
  patterns = get_kmer_abundances_from_kmers.call kmers
  puts [
    seq.definition,
    'start',
    patterns
  ].join("\t")

  # repeat for the end of the contig
  kmers = []
  i = 0
  seq.naseq.reverse_complement.window_search(options[:kmer_size]) do |s|
    kmers.push s.seq.upcase
    i += 1
    break if i>options[:number_of_kmers]
  end
  patterns = get_kmer_abundances_from_kmers.call kmers
  puts [
    seq.definition,
    'end',
    patterns
  ].join("\t")
end

#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'tempfile'
require 'tmpdir'


SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'kmer_abundance_pattern'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <kmer_multiple_abundance_file>

    Given an input kmer then abundances space separated file, and a threshold, print out how many kmers are unique to different subsets of columns\n\n"

  opts.on("--pattern PATTERN", "kmer abundance pattern e.g. '0111001110' [required]") do |arg|
    options[:pattern] = arg.to_i
  end
  opts.on("--kmer-abundances FILE", "kmer multiple abundance file [required]") do |arg|
    options[:abundance_file] = arg
  end
  opts.on("--upper-threshold NUM", "kmer frequency cutoff to saying 'present' [required]") do |arg|
    options[:upper_threshold] = arg.to_i
  end
  opts.on("--lower-threshold NUM", "kmer frequency cutoff to saying 'not present' [required]") do |arg|
    options[:lower_threshold] = arg.to_i
  end
  opts.on("--reads FILES", "comma-separated list of sequence reads files in the same order as the pattern was supplied [required]") do |arg|
    options[:reads_files] = arg.split(',')
  end
  opts.on("--output-dir DIR", "already-created output directory for assembly (currently messy on the output) [required]") do |arg|
    options[:output_directory] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 1 or options[:upper_threshold].nil? or options[:lower_threshold].nil? or options[:pattern].nil? or options[:abundance_file].nil? or options[:reads].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Parse pattern from cmdline
desired_pattern = KmerAbundancePattern.new
desired_pattern.parse_from_human(options[:pattern])
if options[:reads].length != desired_pattern.length
  raise "Number of entries in the pattern #{desired_pattern.length} and number of reads files #{options[:reads].length} not equivalent!"
end

# Collect the kmers that will be used to find trusted reads i.e.
# Go through each line of the kmer abundance file, looking for kmers that suit the pattern
input_file = nil
if options[:kmer_multiple_abundance_file] == '-'
  input_file = $stdin
else
  input_file = options[:kmer_multiple_abundance_file]
end
csv = CSV.new(input_file, :col_sep => ' ')

whitelist_kmers = []
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
      break
    end
  end

  # Reached here means this kmer never fell in no man's land
  if desired_pattern.consistent_with this_pattern
    whitelist_kmers.push row[0]
  end
end
log.info "After parsing the kmer multiple abundance file, found #{whitelist_kmers.length} kmers that matched the pattern"


# grep the pattern out from the raw reads, subsampling so as to not overwhelm the assembler
Tempfile.open('grep_f_input') do |tempfile|
  tempfile.puts whitelist_kmers.join("\n")
  tempfile.close

  threadpool = []
  sampled_read_files = []
  options[:reads_files].each_with_index do |file, i|
    sampled = Tempfile.new('fasta_sample')
    sampled_read_files.push sampled

    next unless desired_pattern[i] #Don't extract reads from reads where those reads should not have been amplified

    #thr = Thread.new do
      grep_cmd = "zcat {file} | fgrep -f #{tempfile.path} -A2 -B1 | grep -v ^--$ | awk '{print \">\" substr($0,2);getline;print;getline;getline}' | fastq_sample -n 200 - > #{sampled.path}"
      log.debug "Running cmd: #{grep_cmd}"
      status, stdout, stderr = systemu grep_cmd
      raise stderr if stderr != ''
      raise unless status.exitstatus == 0
      log.debug "Finished extracting reads from #{file}"
    #end
    threadpool.push thr
  end
  #threadpool.each do |thread| thread.join; end #wait until everything is a-ok


  outdir = options[:output_directory]
  Dir.mkdir outdir unless outdir.exist?
  raise if outdir.directory?

  Dir.chdir outdir do
    log.info "Finished extracting reads for sampling. Now pooling sampled reads"
    pooled_reads_filename = 'pooled_sampled_reads.fasta'
    cat_cmd = "cat #{sampled_read_files.collect{|s| s.path}.join ' '} >#{pooled_reads_filename}"

    log.info "Assembling sampled reads"
    cap3_cmd = "cap3 #{pooled_reads.path}"
    log.debug "Running cmd: #{cap3_cmd}"
    status, stdout, stderr = systemu cap3_cmd
    raise stderr if stderr != ''
    raise unless status.exitstatus == 0
    log.debug "Finished assembly"
  end


  # Clean up
  sampled_read_files.each do |s| s.close; end
end

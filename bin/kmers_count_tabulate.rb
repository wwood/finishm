#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'progressbar'
require 'tempfile'
require 'systemu'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :min_count => 1,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <kmers_count_output1> [<kmers_count_output2> ..]

    Take a list of files output from libngs' kmers_count tool, after being run through gnu sort.

    Create a table, where the columns are each file, the rows are each kmer, and
    the cells are the percent of that file's kmer actually is that kmer.\n\n"


  opts.on("--output-file FILENAME", "Output file path [required]") do |arg|
    options[:output_file] = arg
  end

  opts.on("--percentage", "description [default: #{options[:eg]}]") do
  raise "not yet implemented"
    options[:percentage_outputs] = true
  end
  opts.on("--min-count COUNT", "require at least this many kmers to be output into the output file [default: #{options[:min_count]}]") do |arg|
  raise "not yet implemented"
    options[:min_count] = arg.to_i
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length == 0 or options[:output_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

filenames = ARGV
raise "I need more than 1 file" unless filenames.length > 1
log.info "Joining these files: #{filenames.inspect}"

# run gnu join on each file
current_build_file = filenames[0] #Build off the current build file first, then a tempfile subsequently

Tempfile.open('kmers_join1') do |tempfile1|
  Tempfile.open('kmers_join2') do |tempfile2|
    filenames.each_with_index do |file, i|
      next if i==0

      first_file_output_fields = (2..(i+1)).to_a.collect{|n| "1.#{n.to_s}"}.join(',')
      cmd = "join -a1 -a2 -e 0 -o0,#{first_file_output_fields},2.2 #{current_build_file.inspect} #{file} >#{tempfile2.path}"
      log.info "At #{Time.now}, running #{cmd}.."
      status, stdout, stderr = systemu cmd
      raise stderr unless stderr == ''
      raise 'exitstatus bad1!' unless status.exitstatus == 0
      status, stdout, stderr = systemu "mv #{tempfile2.path} #{tempfile1.path}"
      raise stderr unless stderr == ''
      raise 'exitstatus bad2!' unless status.exitstatus == 0
      current_build_file = tempfile1.path
    end
    status, stdout, stderr = systemu "mv #{current_build_file} #{options[:output_file]}"
    raise stderr unless stderr == ''
    raise 'exitstatus bad3!' unless status.exitstatus == 0
  end
end




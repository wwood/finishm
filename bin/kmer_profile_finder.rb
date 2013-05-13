#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <kmer_multiple_abundance_file>

    Given an input kmer then abundances space separated file, and a threshold, print out how many kmers are unique to different subsets of columns\n\n"

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
if ARGV.length != 1 or options[:upper_threshold].nil? or options[:lower_threshold].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

encoded_counts = {}
max_i =  0

input_file = nil
if ARGV[0] == '-'
  input_file = $stdin
else
  input_file = File.open ARGV[0]
end
csv = CSV.new(input_file, :col_sep => ' ')

csv.each do |row|
  kmer = row[0]
  counts = row[1...row.length].collect{|s| s.to_i}
  index = 0
  counts.each_with_index do |count, i|
    max_i = i if i > max_i

    if count > options[:upper_threshold]
      increment = (1<<i)
      index += increment
      log.debug "Found a passable for #{options[:threshold]} in index #{i} for #{counts}, count is now #{index}" if log.debug?
    elsif count < options[:lower_threshold]
      # do nothing
    else
      # coverage was in no man's land between thresholds.
      # Ignore this kmer as noise.
      break
    end
  end

  if index != 0
    encoded_counts[index] ||= 0
    encoded_counts[index] += 1
  end
end

(0..encoded_counts.keys.max).each do |i|
  total = encoded_counts[i]
  unless total.nil?
    unencoded = i.to_s(2)

    while unencoded.length <= max_i
      unencoded = '0'+unencoded
    end

    puts [
      i,
      total,
      unencoded,
    ].join "\t"
  end
end

#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :min => 0,
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    grep a multiple kmer abundance file according to specified criteria\n\n"

  opts.on("--min NUMBER", "At least 1 column has at least this many observations [default: #{options[:min]}]") do |arg|
    options[:min] = arg.to_f
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 1
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

CSV.foreach(ARGV[0], :col_sep => ' ') do |row|
  kmer = row[0]
  passable = false
  row[1...row.length].each do |count|
    if count.to_f > options[:min]
      passable = true
      break
    end
  end
  puts row.join(' ') if passable
end


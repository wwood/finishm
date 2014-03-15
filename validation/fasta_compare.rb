#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Are the sequences in two files the same?\n\n"

  opts.on("-1", "--fasta-1 FILE", "path to the first fasta file [required]") do |arg|
    options[:first] = arg
  end
  opts.on("-2", "--fasta-2 FILE", "path to the second fasta file [required]") do |arg|
    options[:second] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:first].nil? or options[:second].nil?
  $stderr.puts "More arguments required"
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')


#hash the first file
firsts = {}
Bio::FlatFile.foreach(options[:first]) do |seq|
  name = seq.definition
  raise "Duplicated seq name #{name}" if firsts.key?(name)

  firsts[name] = seq.seq
end

sames = []
differents = []
Bio::FlatFile.foreach(options[:second]) do |seq|
  name = seq.definition

  unless firsts.key?(name)
    raise "Unable to find sequence #{name} in the first file, when it exists in the second"
  end
  if firsts[name] == seq.seq
    sames.push name
  else
    differents.push name
  end
end
log.info "Found #{sames.length} sequences that were the same and #{differents.length} that were different"
if log.debug?
  sames.each do |same|
    log.debug "#{same} was the same in each"
  end
  differents.each do |d|
    log.debug "#{d} was different"
  end
end

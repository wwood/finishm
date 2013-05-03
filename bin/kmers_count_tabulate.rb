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
    Usage: #{SCRIPT_NAME} <kmers_count_output1> [<kmers_count_output2> ..]

    Take a list of files output from libngs' kmers_count tool, after being run through gnu sort.

    Create a table, where the columns are each file, the rows are each kmer, and
    the cells are the percent of that file's kmer actually is that kmer.\n\n"


  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length == 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


files = ARGV.collect do |arg|
  CSV.open arg, :col_sep => ' '
end

current_rows = files.collect do |f|
  r = f.shift
end

# Print headers
print "\t"
puts ARGV.collect{|arg| File.basename arg}.join("\t")

all_finished = false
finished = [false]*files.length

while !all_finished
  # Choose the lowest lexigraphical kmer from the files
  lowest_kmer = current_rows.collect { |r|
    if r.nil?
      'Z'
    else
      r[0]
    end
  }.min

  print lowest_kmer
  files.each_with_index do |f, i|
    print "\t"
    if !finished[i] and current_rows[i][0] == lowest_kmer
      print current_rows[i][1]
      current_rows[i] = files[i].shift
      if current_rows[i].nil?
        finished[i] = true
      else
        raise "Unexpected input format with this line: #{current_rows[i]}" unless current_rows[i].length == 2
      end
    else
      print 0
    end
  end
  puts

  all_finished = finished.reject{|f| f}.length == 0
end

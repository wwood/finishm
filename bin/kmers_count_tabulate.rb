#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'progressbar'

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

  opts.on("--percentage", "description [default: #{options[:eg]}]") do
    options[:percentage_outputs] = true
  end
  opts.on("--min-count COUNT", "require at least this many kmers to be output into the output file [default: #{options[:min_count]}]") do |arg|
    options[:min_count] = arg.to_i
  end

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

filenames = ARGV
files = filenames.collect do |arg|
  CSV.open arg, :col_sep => ' '
end

current_rows = files.collect do |f|
  r = f.shift
end

total_counts = []
line_count0 = 0
progress = nil
if options[:percentage_outputs]
  log.info "First pass - reading total number of kmers in each sample"
  first_file = true
  total_counts = filenames.collect do |f|
    count = 0
    CSV.foreach(f, :col_sep => ' ') do |row|
      count += row[1].to_i
      line_count0 += 1 if first_file
    end
    first_file = false
    count
  end
  log.info "Found #{filenames.length} files to kmer count, and there's #{line_count0} lines in the first one"
  log.info "Found total different kmers in each one: #{total_counts.inspect}"
  progress = ProgressBar.new('kmer_counting',line_count0) if log.info?
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

  total_observations_of_this_kmer = 0
  to_print = lowest_kmer
  files.each_with_index do |f, i|
    to_print += "\t"
    if !finished[i] and current_rows[i][0] == lowest_kmer
      total_observations_of_this_kmer += current_rows[i][1].to_i
      if options[:percentage_outputs]
        to_print += (current_rows[i][1].to_f / total_counts[i]).to_s
      else
        to_print += current_rows[i][1]
      end
      current_rows[i] = files[i].shift
      if current_rows[i].nil?
        finished[i] = true
      else
        raise "Unexpected input format with this line: #{current_rows[i]}" unless current_rows[i].length == 2
        unless progress.nil?
          progress.inc if i==0 and current_rows[i][0] == lowest_kmer
        end
      end
    else
      to_print += '0'
    end
  end

  if !options[:min_count] or total_observations_of_this_kmer >= options[:min_count]
    puts to_print
  end

  all_finished = finished.reject{|f| f}.length == 0
end
progress.finish unless progress.nil?

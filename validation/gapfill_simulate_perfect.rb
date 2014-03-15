#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'bio-commandeer'
require 'set'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :max_gap_size => 1000,
  :min_gap_size => 100,
  :contig_end_buffer => 2000,
  :num_to_gap => 100,

  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

Introduce random gaps into sequences for the purpose of validating gap filling algorithms.\n\n"

  opts.on("-i", "--input-fasta FILE", String, "Input fasta file of sequences for simulation [required]") do |arg|
    options[:input_fasta_file] = arg
  end
  opts.on("--gapped FILE", String, "Output file with simulated gaps [required]") do |arg|
    options[:gapped_output_file] = arg
  end
  opts.on("--answer FILE", String, "Output file with correct sequences, not including sequences that did not have gaps introduced [required]") do |arg|
    options[:answer_output_file] = arg
  end

  opts.separator "\nOptional parameters:\n\n"
  opts.on("--num-to-gap NUM", Integer, "Number of sequences to introduce gaps into [default: #{options[:num_to_gap]}]") do |arg|
    options[:num_to_gap] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME); log.outputters[0].formatter = Log4r::PatternFormatter.new(:pattern => "%5l %c %d: %m", :date_pattern => '%d/%m %T')

# Assumes that files to be simulated are >10kb in length
# Count the number of sequences in the file with grep
# choose 100 random ones to have gaps introduced artificially
cmd = "cat #{options[:input_fasta_file]} |grep -c '>'"
num_seqs = Bio::Commandeer.run cmd, :log => log
num_seqs = num_seqs.to_i
log.info "Found #{num_seqs} sequences in input file"
raise "Input file problem - not enough seqs" if num_seqs < options[:num_to_gap]

# Choose 100 random sequences, choosing indices
indices_to_gap = Set.new ((0...num_seqs).to_a.sample(options[:num_to_gap]))

# Read in fasta file into big hash, keeping only those sequences that have
# been chosen randomly.
i=0
randomer = Random.new
gapped_output = File.open(options[:gapped_output_file],'w')
answer_output = File.open(options[:answer_output_file],'w')
j=0
Bio::FlatFile.foreach(options[:input_fasta_file]) do |seq|
  if indices_to_gap.include?(i)
    # introduce 100-1000bp gaps - choose length from uniform distribution
    gap_size = randomer.rand(options[:max_gap_size]-options[:min_gap_size] )+options[:min_gap_size]
    gap_size = gap_size.to_i

    # ensure that the gap is more than 1kb away from the end of the sequence
    seq_length = seq.seq.length
    length_for_simulation = seq_length-2*options[:contig_end_buffer]-gap_size
    raise "sequence length problem in #{seq.definition} - not long enough" unless length_for_simulation > 0
    offset = randomer.rand(length_for_simulation).to_i

    # log
    name = ">#{j} #{seq.definition}"
    log.info "Introducing a gap fo size #{gap_size} at position #{offset} into #{name}"

    # Print out the gapped with introduced gaps
    gapped_output.puts name
    gapped_output.print seq.seq[0...offset]
    gapped_output.print 'N'*gap_size
    gapped_output.puts seq.seq[(offset+gap_size)..-1]

    # Print out the answer without introduced gaps
    answer_output.puts name
    answer_output.puts seq.seq

    j += 1
  end
  i += 1
end

gapped_output.close
answer_output.close



#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio-samtools'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :contig_end_length => 200,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>

      Takes a sorted, indexed BAM file and outputs the coverages of each of the reference sequences\n\n"

    opts.on("-b", "--bam BAM_FILE", "BAM file that defines overall mapping/coverage [required]") do |arg|
      options[:bam_file] = arg
    end
    opts.on("-f", "--reference-fasta FASTA_FILE", "FASTA file of the reference [required]") do |arg|
      options[:fasta_file] = arg
    end

    opts.on("-l", "--end-length LENGTH", "How far from the end to count [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
      raise "Inappropriate end length detected, I need a positive number, found #{options[:contig_end_length]} parsed from #{arg}" if options[:contig_end_length] < 1
    end

    # logger options
    opts.separator "\n\tVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0 or options[:bam_file].nil? or options[:fasta_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


  # open the BAM file for reading
  bam = Bio::DB::Sam.new(:bam => options[:bam_file], :fasta => options[:fasta_file])
  bam.open
  puts %w(Reference StartCoverage EndCoverage).join("\t")
  contigs_file = Bio::FlatFile.auto(options[:fasta_file])
  contigs_file.each_entry do |record|
    ref = record.definition
    ref_length = record.length
  #bam.each_reference do |ref, ref_length| #currently commented out because there is extraneous output on STDOUT if you use this
    next if ref == '*' #Ignore umapped reads
    # Coverage of the start
    end_length = options[:contig_end_length]
    end_length = ref_length if ref_length < options[:contig_end_length]

    cov_start = bam.average_coverage(ref, 1, end_length)
    cov_end = bam.average_coverage(ref, ref_length-end_length+1, end_length)
    puts [
      ref, cov_start, cov_end
    ].join("\t")
  end
  exit
end #end if running as a script

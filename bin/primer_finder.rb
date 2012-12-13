#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'progressbar'

$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'priner'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :reverse_complement => false,
    :gc_clamp_length => 2,
    :min_melting_temperature => 56,
    :max_melting_temperature => 62,
    :min_primer_length => 15,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments> <sequence>
      
      Take a sequence, and find an oligos that fit certain parameters. Sort of like primer3_core but only look for single oligos, not pairs.\n\n"

    opts.on("-s", "--sequence-file FILE", "Fasta file of sequences [required]") do |arg|
      options[:input_file] = arg 
    end

    opts.separator "\nOptional arguments:\n\n"  
    opts.on("-r", "--reverse-complement", "Design primers pointing backwards off the start of the sequence in reverse, not off the end forwards [default: #{options[:reverse_complement]}]") do
      options[:reverse_complement] = true
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  # Read the contigs in
  contigs = []
  Bio::FlatFile.foreach(options[:input_file]) do |entry|
    contigs.push entry.seq.to_s
  end
  log.info "Read in #{contigs.length} contigs from #{options[:contigs_file]}"
  raise unless contigs.length == 1
  
  seq = contigs[0]
  
  # Reverse complement if required
  if options[:reverse_complement]
    seq = Bio::Sequence::NA.new(seq).reverse_complement.to_s
  end
  seq.upcase!
  
  raise "Whackiness in the supplied sequence! #{seq}" unless seq.match(/^[ATGC]+$/)
  raise unless options[:min_primer_length] > options[:gc_clamp_length]
  
  # Find out all those positions that have a GC clamp of enough
  gc_clamped_positions = []
  (0...seq.length).each do |pos|
    next unless pos-options[:min_primer_length] >= -1
    
    if seq[pos-options[:gc_clamp_length]+1..pos].match(/^[GC]+$/)
      gc_clamped_positions << pos
    end
  end
  log.info "Found #{gc_clamped_positions.length} positions with a suitable GC-clamp"
  
  # Find those with suitable melting temperatures
  designer = OligoDesigner.new
  progress = ProgressBar.new('primer_finding', gc_clamped_positions.length)
  gc_clamped_positions.each do |pos|
    # Iteratively make primers longer. Start with min primer length, end when the Tm exceeds the maximum allowable
    current_length = options[:min_primer_length]
    
    over_tm_max = false
    while !over_tm_max and pos-current_length >= 0
      oligo = seq[pos-current_length+1..pos]
      tm = designer.melting_temperature oligo
      puts oligo
      puts tm
      exit
      
      if tm > options[:min_melting_temperature]
        if tm < options[:max_melting_temperature]
          # This is a hit
          cur_seq = oligo
          if options[:reverse_complement]
            cur_seq = Bio::Sequence::NA.new(cur_seq).reverse_complement.to_s
          end
          
          puts [
            cur_seq, tm
          ].join("\t")
        else
          over_tm_max =  true
          break #making the primer longer can only result in higher melting temperatures, and we are already over the max
        end
      end
      current_length += 1
    end
    progress.inc
  end
  progress.finish
  
  
end #end if running as a script
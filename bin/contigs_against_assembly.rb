#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'systemu'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} -b <contigs_against_assembly.blast_outfmt6.csv>

    Takes a set of contigs, and an assembly. Works out if there are any contigs where there is a blast hit spanning of the contigs using two of the assembly's contig ends.\n\n"

  opts.on("--query FASTA_FILE", "new contigs fasta file [Required]") do |arg|
    options[:query_file] = arg
  end
  opts.on("--blastdb FASTA_FILE_FORMATTED", "basename of makeblastdb output [Required]") do |arg|
    options[:blastdb] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:query_file].nil? or options[:blastdb].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Read in the blast file
blast_results = []
class BlastResult
  attr_accessor :qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :subject_end, :evalue, :bitscore, :query_length, :subject_length

  attr_accessor :cutoff_inwards

  def initialize
    @cutoff_inwards = 500
  end

  def hits_end_of_subject?
    @subject_end >= @subject_length-@cutoff_inwards and @length >= 100
  end

  def hits_start_of_subject?
    @sstart <= @cutoff_inwards and @length >= 100
  end

  def hits_end_of_query?
    @qend >= @query_length-@cutoff_inwards and @length >= 100
  end

  def hits_start_of_query?
    @qstart <= @cutoff_inwards and @length >= 100
  end
end

status, blast_output, stderr = systemu "blastn -query #{options[:query_file].inspect} -db #{options[:blastdb].inspect} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue 1e-5"
raise stderr unless stderr==""
raise "bad status running blast" unless status.exitstatus == 0
log.debug "Finished running blast, presumably successfully"

blast_output.each_line do |line|
  res = BlastResult.new
  row = line.chomp.split "\t"
  [:qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart,
  :qend, :sstart, :subject_end, :evalue, :bitscore,
  :query_length, :subject_length].each_with_index do |attr, i|
    res.send "#{attr}=".to_sym, row[i]
  end
  [:length, :mismatch, :gapopen, :qstart,
  :qend, :sstart, :subject_end,:query_length, :subject_length].each do |attr|
    res.send "#{attr}=".to_sym, res.send(attr).to_i
  end
  [:pident, :evalue, :bitscore].each do |attr|
    res.send "#{attr}=".to_sym, res.send(attr).to_f
  end

  blast_results.push res
end
log.info "Parsed #{blast_results.length} blast results e.g. #{blast_results[0].inspect}"


query_to_blast_results = {}
hit_to_blast_results = {}
blast_results.each do |result|
  query_to_blast_results[result.qseqid] ||= []
  query_to_blast_results[result.qseqid].push result

  hit_to_blast_results[result.sseqid] ||= []
  hit_to_blast_results[result.sseqid].push result
end

# For each query sequence, does it map to the ends of both contigs
header = %w(query subject1 subject2 qstart1? qend1? sstart1? send1? qstart2? qend2? sstart2? send2?).join("\t")
query_to_blast_results.each do |query_id, hits|
  query_length = hits[0].query_length
  keepers = []

  hits.each do |hit|
    # perfect if it hits the start or the end (but not both) of both the query and the subject, unless it is circular
    if hit.hits_start_of_query? ^ hit.hits_end_of_query? and
      hit.hits_start_of_subject? ^ hit.hits_end_of_subject?
      keepers.push hit
    elsif hit.hits_start_of_query? or hit.hits_end_of_query? or
      hit.hits_start_of_subject? or hit.hits_end_of_subject?
      log.info "There's a half-correct hit for #{query_id}: qstart? #{hit.hits_start_of_query?} qend #{hit.hits_end_of_query?} "+
        "sstart #{hit.hits_start_of_subject?} send #{hit.hits_end_of_subject?}, to subject sequence #{hit.sseqid}"
    end
  end

  if keepers.empty?
    log.debug "no latchings found for #{query_id}"
  elsif keepers.length == 1
    log.info "Query #{query_id} only latches on to a single end, maybe manually inspect"
  elsif keepers.length == 2
    log.debug "Query #{query_id} has 2 keepers!"
    q = keepers.collect{|hit| hit.hits_start_of_query?}.join
    s = keepers.collect{|hit| hit.hits_start_of_subject?}.join
    if (q == 'truefalse' or q == 'falsetrue') and
      (s == 'truefalse' or s == 'falsetrue')
      outs = (0..1).collect{|i|
        [
          keepers[i].hits_start_of_query?,
          keepers[i].hits_end_of_query?,
          keepers[i].hits_start_of_subject?,
          keepers[i].hits_end_of_subject?,
        ]
      }.flatten
      unless header.nil?
        puts header
        header = nil
      end
      puts [query_id, keepers[0].sseqid, keepers[1].sseqid, outs].flatten.join("\t")
    else
      log.info "Query #{query_id} has 2 keepers, but they are fighting it seems"
    end
  else
    log.info "More than 2 keepers found for #{query_id}, manual inspection likely required"
  end
end


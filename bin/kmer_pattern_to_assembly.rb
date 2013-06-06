#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'tempfile'
require 'pp'
require 'systemu'
require 'bio-velvet'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'finishm'
$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'priner'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
  :min_leftover_length => false,
  :kmer_coverage_target => 1,
  :velvet_kmer_size => 155,
  :terminal_contig_search_kmer_size => 33,
  :contig_end_length => 300,
  :graph_search_leash_length => 20000,
  :reads_to_assemble => nil,
}

# TODO: make a better interface for this. Maybe specify an entire genome, and then "Contig_1 end, Contig_3 start" or something
# Look at the last 300bp of the first contig.
extract_exactly_one_contig_from_file = lambda do |fasta_file_path|
  contig = nil
  Bio::FlatFile.foreach(Bio::FastaFormat, fasta_file_path) do |e|
    if contig.nil?
      contig = e.seq
    else
      raise "Multiple sequences found in a contig file! I need exactly one"
    end
  end
  raise "I need a contig to be in the start contig file" if contig.nil?
  Bio::Sequence::NA.new(contig.to_s)
end

o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <kmer_multiple_abundance_file>

    Given an input kmer then abundances space separated file, and a threshold, print out how many kmers are unique to different subsets of columns\n\n"

  opts.on("--pattern PATTERN", "kmer abundance pattern e.g. '0111001110' [required]") do |arg|
    options[:pattern] = arg
  end
  opts.on("--kmer-abundances FILE", "kmer multiple abundance file [required]") do |arg|
    options[:kmer_multiple_abundance_file] = arg
  end
  opts.on("--upper-threshold NUM", "kmer frequency cutoff to saying 'present' [required]") do |arg|
    options[:upper_threshold] = arg.to_i
  end
  opts.on("--lower-threshold NUM", "kmer frequency cutoff to saying 'not present' [required]") do |arg|
    options[:lower_threshold] = arg.to_i
  end
  opts.on("--reads FILES", "comma-separated list of sequence reads files in the same order as the pattern was supplied [required]") do |arg|
    options[:reads_files] = arg.split(',').collect{|r| File.absolute_path r}
  end
  opts.on("--start-contig FASTA", "path to a fasta file with the starting contig in it (only). Assumes we are building off the end of this contig [required]") do |arg|
    options[:start_contig] = extract_exactly_one_contig_from_file.call arg
  end
  opts.on("--end-contig FASTA", "path to a fasta file with the ending contig in it (only). Assumes we are building onto the start of this contig [required]") do |arg|
    options[:end_contig] = extract_exactly_one_contig_from_file.call arg
  end

  opts.separator "\nOptional arguments:\n\n"
  opts.on("--min-leftover-read-length NUMBER", "when searching for reads with kmers, require the kmer to be at the beginning or end of the selected read [default: #{options[:min_leftover_length]}]") do |arg|
    options[:min_leftover_length] = arg.to_i
  end
  opts.on("--kmer-coverage-target NUMBER", "when searching for reads with kmers, require this many copies per kmer [default: #{options[:kmer_coverage_target]}]") do |arg|
    options[:kmer_coverage_target] = arg.to_i
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:upper_threshold].nil? or options[:lower_threshold].nil? or options[:pattern].nil? or options[:kmer_multiple_abundance_file].nil? or options[:reads_files].nil?
pp options
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
Bio::Log::LoggerPlus.new 'bio-velvet'
Bio::Log::CLI.configure 'bio-velvet'

if(false)
# Parse pattern from cmdline
desired_pattern = KmerAbundancePattern.new
desired_pattern.parse_from_human(options[:pattern])
if options[:reads_files].length != desired_pattern.length
  raise "Number of entries in the pattern #{desired_pattern.length} and number of reads files #{options[:reads].length} not equivalent!"
end

# Collect the kmers that will be used to find trusted reads i.e.
# Go through each line of the kmer abundance file, looking for kmers that suit the pattern
input_file = nil
if options[:kmer_multiple_abundance_file] == '-'
  input_file = $stdin
else
  input_file = File.open options[:kmer_multiple_abundance_file]
end
csv = CSV.new(input_file, :col_sep => ' ')

whitelist_kmers = []
blacklist_kmers = []
csv.each do |row|
  max_i = row.length - 2 if max_i.nil?

  kmer = row[0]
  counts = row[1...row.length].collect{|s| s.to_i}

  this_pattern = []
  counts.each_with_index do |count, i|
    if count > options[:upper_threshold]
      this_pattern[i] = true
    elsif count < options[:lower_threshold]
      this_pattern[i] = false
    else
      # coverage was in no man's land between thresholds.
      # Ignore this kmer as noise.
      break
    end
  end
  #log.debug "Found pattern #{this_pattern} from kmer #{kmer}, which has abundances #{counts}" if log.debug?
  next unless this_pattern.length == desired_pattern.length

  # Reached here means this kmer never fell in no man's land
  if desired_pattern.consistent_with? this_pattern
    whitelist_kmers.push row[0]
  else
    # kmer is not present when it should be
    blacklist_kmers.push row[0]
  end
end
log.info "After parsing the kmer multiple abundance file, found #{whitelist_kmers.length} kmers that matched the pattern, and #{blacklist_kmers.length} that didn't"
unless whitelist_kmers.length > 0
  log.error "No kmers found that satisfy the given pattern, exiting.."
  exit 1
end


#outdir = options[:output_directory]
#Dir.mkdir outdir unless Dir.exist?(outdir)

# grep the pattern out from the raw reads, subsampling so as to not overwhelm the assembler
pooled_reads_filename = 'pooled_sampled_reads.fasta'
#Tempfile.open('whitelist') do |white|
File.open 'whitelist', 'w' do |white|
  white.puts whitelist_kmers.join("\n")
  white.close

  #Tempfile.open('blacklist') do |black|
  File.open('black','w') do |black|
    black.puts blacklist_kmers.join("\n")
    black.close

    threadpool = []
    sampled_read_files = []
    options[:reads_files].each_with_index do |file, i|
      next unless desired_pattern[i] #Don't extract reads from reads where those reads should not have been amplified

      sampled = File.basename(file)+'.sampled_reads.fasta'
      sampled_read_files.push sampled

      grep_path = "#{ENV['HOME']}/git/priner/bin/read_selection_by_kmer "
      if options[:min_leftover_length]
        grep_path += "--min-leftover-length #{options[:min_leftover_length]} "
      end
      thr = Thread.new do
        grep_cmd = "#{grep_path} --whitelist #{white.path} --blacklist #{black.path} --reads #{file} --kmer-coverage-target #{options[:kmer_coverage_target]} > #{sampled}"
        log.debug "Running cmd: #{grep_cmd}"
        status, stdout, stderr = systemu grep_cmd
        log.debug stderr

        raise unless status.exitstatus == 0
        log.debug "Finished extracting reads from #{file}"
      end
      threadpool.push thr
    end
    threadpool.each do |thread| thread.join; end #wait until everything is finito

    log.info "Finished extracting reads for sampling. Now pooling sampled reads"
    pool_cmd = "cat #{sampled_read_files.join ' '} >#{pooled_reads_filename}"
    log.debug "Running cmd: #{pool_cmd}"
    status, stdout, stderr = systemu pool_cmd
    raise stderr if stderr != ''
    raise unless status.exitstatus == 0
  end
end

end

pooled_reads_filename = 'pooled_sampled_reads.fasta'
log.info "Assembling sampled reads with velvet"
velvet_result = Bio::Velvet::Runner.new.velvet(options[:velvet_kmer_size], "-short #{pooled_reads_filename}", '-cov_cutoff 1.5')
log.info "Finished running assembly"

log.info "Parsing the graph output from velvet"
graph = velvet_result.last_graph
pp graph.arcs
log.info "Finished parsing graph, and found #{graph.nodes.length} nodes"

log.info "Finding kmers that are specific to the end of the first contig"
start_contig = options[:start_contig]
end_contig = options[:end_contig]
if [start_contig.length, end_contig.length].min < 2*options[:contig_end_length]
  # TODO: if the contig is very short, earlier kmers in the kmers array may
  # not be closer to the middle of the contig. re-order the input kmer hash to make it so?
  log.warn "Choice of initial/terminal nodes to perform graph search with may not be optimal due to the small contig size"
end

start_kmers = []
len = start_contig.length
start_contig.subseq(len-options[:contig_end_length],len).window_search(options[:terminal_contig_search_kmer_size]) do |kmer_na|
  start_kmers.push kmer_na.to_s.upcase
end
end_kmers = []
end_contig.subseq(1,options[:contig_end_length]).reverse_complement.window_search(options[:terminal_contig_search_kmer_size]) do |kmer_na|
  end_kmers.push kmer_na.to_s.upcase
end

finder = Bio::AssemblyGraphAlgorithms::NodeFinder.new
start_node, start_node_forward = finder.find_unique_node_with_kmers(graph, start_kmers)
log.info "Finding kmers that are specific to the start of the second contig"
end_node, end_node_forward = finder.find_unique_node_with_kmers(graph, end_kmers)
if start_node.nil? or end_node.nil?
  log.error "Unable to find any nodes in the graph that have suitable kmers in them, sorry. Maybe fix the node finding code?"
  exit
end
log.info "Node(s) found that are suitable as initial and terminal nodes in the graph search, respectively: #{start_node.node_id} and #{end_node.node_id}"

log.info "Searching for trails between the initial and terminal nodes, within the assembly graph"
cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
trails = cartographer.find_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
log.info "Found #{trails.length} trail(s) between the initial and terminal nodes"

log.debug "Found trails: #{trails.collect{|t| "Trail: #{t.collect{|n| n.node_id}.join(',')}"}.join(', ')}"

sequencer = Bio::AssemblyGraphAlgorithms::LazyGraphWalker.new
trails.each_with_index do |trail,i|
  seq = sequencer.trail_sequence(graph,trail)
  puts ">finishm_trail_#{i+1}"
  puts seq
end

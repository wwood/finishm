#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'tempfile'
require 'pp'
require 'systemu'
require 'bio-velvet'
require 'set'

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
  :contig_end_length => 300,
  :graph_search_leash_length => 20000,
  :reads_to_assemble => nil,
  :assembly_coverage_cutoff => 1.5,
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
  opts.on("--assembly-kmer NUMBER", "when assembling, use this kmer length [default: #{options[:velvet_kmer_size]}]") do |arg|
    options[:velvet_kmer_size] = arg.to_i
  end
  opts.on("--already-patterned-reads FILE", "Attempt to assemble the reads in the specified file, useful for re-assembly [default: off]") do |arg|
    options[:already_patterned_reads] = arg
  end
  opts.on("--output-assembly PATH", "Output assembly intermediate files to this directory [default: off]") do |arg|
    options[:output_assembly_path] = arg
  end
  opts.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
    options[:output_graph_png] = arg
  end
  opts.on("--assembly-svg PATH", "Output assembly as an SVG file [default: off]") do |arg|
    options[:output_graph_svg] = arg
  end
  opts.on("--assembly-dot PATH", "Output assembly as an DOT file [default: off]") do |arg|
    options[:output_graph_dot] = arg
  end
#  opts.on("--output-begin-kmers PATH", "Output kmers found at the beginning point to this file [default: off]") do |arg|
#    options[:output_begin_kmers] = arg
#  end
#  opts.on("--output-end-kmers PATH", "Output kmers found at the ending point to this file [default: off]") do |arg|
#    options[:output_end_kmers] = arg
#  end
  opts.on("--assembly-coverage-cutoff NUMBER", "Require this much coverage in each node, all other nodes are removed [default: #{options[:assembly_coverage_cutoff]}]") do |arg|
    options[:assembly_coverage_cutoff] = arg.to_f
  end
  opts.on("--contig-end-length LENGTH", "Number of base pairs to start into the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
    options[:contig_end_length] = arg.to_i
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

pooled_reads_filename = 'pooled_sampled_reads.fasta'
if options[:already_patterned_reads] #If skipping read extraction
  pooled_reads_filename = options[:already_patterned_reads]

else
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
probe = 'TTACATCTTATCTACAATAAACCTTCTGCCTTAGTTTTAGAGCCTATCCGAAAAGTCCTGCTGCTCTGAATGTTATCCAAGCACATGCAAAATGAATTAGT'
    this_pattern = []
    counts.each_with_index do |count, i|
      if count > options[:upper_threshold]
        this_pattern[i] = true
      elsif count < options[:lower_threshold]
        this_pattern[i] = false
      else
        # coverage was in no man's land between thresholds.
        # Ignore this kmer as noise.
        this_pattern[i] = '-'
      end
    end
    #log.debug "Found pattern #{this_pattern} from kmer #{kmer}, which has abundances #{counts}" if log.debug?

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
      log.info "Extracting reads that contain suitable kmers"
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

log.info "Extracting dummy reads from the ends of contigs to use as anchors"
start_contig = options[:start_contig]
end_contig = options[:end_contig]
if [start_contig.length, end_contig.length].min < 2*options[:contig_end_length]
  log.warn "Choice of initial/terminal nodes to perform graph search with may not be optimal due to the small contig size"
end
if [start_contig.length, end_contig.length].min < options[:contig_end_length]
  log.error "At least one contig too small to proceed with current code base, need to fix the code to allow such a small contig"
  exit 1
end
# Use the last bit of the first contig and the first bit of the second contig as the anchors
velvet_result = nil
Tempfile.open('anchors.fa') do |tempfile|
  # Putting these same sequences in many times seems to better the
  # chances velvet won't throw them out
  50.times do
    tempfile.puts ">start_contig"
    tempfile.puts start_contig[start_contig.length-options[:contig_end_length]...start_contig.length]
    tempfile.puts ">end_contig"
    #Have to be in reverse, because the node finder finds the node at the start of the read, not the end
    fwd2 = Bio::Sequence::NA.new(end_contig[0...options[:contig_end_length]])
    tempfile.puts fwd2.reverse_complement.to_s.gsub(/n+$/,'')
  end
  tempfile.close
  puts `cat #{tempfile.path}`

  log.info "Assembling sampled reads with velvet"
  # Bit of a hack, but have to use -short1 as the anchors because then start and end anchors will have node IDs 1 and 2, respectively.
  velvet_result = Bio::Velvet::Runner.new.velvet(
    options[:velvet_kmer_size],
    "-short #{tempfile.path} -short2 #{pooled_reads_filename}",
    "-cov_cutoff 0 -read_trkg yes",
    :output_assembly_path => options[:output_assembly_path]
  )
  if log.debug?
    log.debug "velveth stdout: #{velvet_result.velveth_stdout}"
    log.debug "velveth stderr: #{velvet_result.velveth_stderr}"
    log.debug "velvetg stdout: #{velvet_result.velvetg_stdout}"
    log.debug "velvetg stderr: #{velvet_result.velvetg_stderr}"
  end
  log.info "Finished running assembly"
end

log.info "Parsing the graph output from velvet"
graph = Bio::Velvet::Graph.parse_from_file(File.join velvet_result.result_directory, 'Graph2')
log.info "Finished parsing graph: found #{graph.nodes.length} nodes and #{graph.arcs.length} arcs"

if options[:assembly_coverage_cutoff]
  log.info "Removing low-coverage nodes from the graph (less than #{options[:assembly_coverage_cutoff]})"
  cutoffer = Bio::AssemblyGraphAlgorithms::CoverageBasedGraphFilter.new
  deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, options[:assembly_coverage_cutoff], :whitelisted_sequences => [1,2])

  log.info "Removed #{deleted_nodes.length} nodes and #{deleted_arcs.length} arcs from the graph due to low coverage"
  log.info "Now there is #{graph.nodes.length} nodes and #{graph.arcs.length} arcs remaining"
end

finder = Bio::AssemblyGraphAlgorithms::NodeFinder.new
log.info "Finding node representing the end of the first contig"
start_node, start_node_forward = finder.find_unique_node_with_sequence_id(graph, 1)
log.info "Finding node representing the start of the second contig"
end_node, end_node_forward = finder.find_unique_node_with_sequence_id(graph, 2)#TODO: find the node nearest the end of this, not the start
if start_node.nil? or end_node.nil?
  if start_node.nil?
    log.error "Unable to find any nodes in the graph that have kmers corresponding to the _start_ point in them, sorry. Maybe fix the node finding code?"
  end
  if end_node.nil?
    log.error "Unable to find any nodes in the graph that have kmers corresponding to the _end_ point in them, sorry. Maybe fix the node finding code?"
  end

  if options[:output_graph_png] or options[:output_graph_svg] or options[:output_graph_dot]
    log.info "Converting assembly to a graphviz PNG/SVG/DOT, even if start/end node was not be found properly"
    viser = Bio::Assembly::ABVisualiser.new
    gv = viser.graphviz(graph)
    if options[:output_graph_png]
      log.info "Writing PNG of graph to #{options[:output_graph_png]}"
      gv.output :png => options[:output_graph_png]
    end
    if options[:output_graph_svg]
      log.info "Writing SVG of graph to #{options[:output_graph_svg]}"
      gv.output :svg => options[:output_graph_svg]
    end
    if options[:output_graph_dot]
      log.info "Writing DOT of graph to #{options[:output_graph_dot]}"
      gv.output :dot => options[:output_graph_dot]
    end
  end
  log.error "Unknown start or end points, giving up, sorry."
  exit 1
end
log.info "Node(s) found that are suitable as initial and terminal nodes in the graph search, respectively: #{start_node.node_id} and #{end_node.node_id}"

if options[:output_graph_png]
  log.info "Converting assembly to a graphviz PNG"
  viser = Bio::Assembly::ABVisualiser.new
  gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id})
  gv.output :png => options[:output_graph_png], :use => :neato
end
if options[:output_graph_svg]
  log.info "Converting assembly to a graphviz SVG"
  viser = Bio::Assembly::ABVisualiser.new
  gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id})
  gv.output :svg => options[:output_graph_svg], :use => :neato
end
if options[:output_graph_dot]
  log.info "Converting assembly to a graphviz DOT"
  viser = Bio::Assembly::ABVisualiser.new
  gv = viser.graphviz(graph, {:start_node_id => start_node.node_id, :end_node_id => end_node.node_id, :digraph => false})
  gv.output :dot => options[:output_graph_dot]
end

log.info "Searching for trails between the initial and terminal nodes, within the assembly graph"
cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#raise "Untested connection finder below"
#trails = cartographer.find_all_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
trails = cartographer.find_trails_between_nodes(graph, start_node, end_node, options[:graph_search_leash_length], start_node_forward)
log.info "Found #{trails.length} trail(s) between the initial and terminal nodes"

log.debug "Found trails: #{trails.collect{|t| t.to_s}.join("\n")}"

trails.each_with_index do |trail, i|
  puts ">trail#{i+1}"
  puts trail.sequence
end

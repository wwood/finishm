#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'pp'
require 'bio'
require 'pry'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

$:.unshift File.join(ENV['HOME'],'git/yargraph/lib')
require 'yargraph'

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Description of what this program does...\n\n"

  opts.on("-c", "--connnections FILE", "connections file [required]") do |arg|
    options[:connections_file] = arg
  end
  opts.on("-f", "--fasta FILE", "Fasta file of all contigs [required]") do |arg|
    options[:fasta_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:connections_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)

class Probe
  attr_accessor :contig_name, :side

  def to_setable
    [@contig_name, @side]
  end
end

class ProbeSet

end

graph = Yargraph::UndirectedGraph.new

num_probes_circular = 0
CSV.foreach(options[:connections_file], :col_sep => "\t") do |row|
  if row.length == 3
    # e.g. 110811_E_1_D_nesoni_single_contig_1011407.start	110811_E_1_D_nesoni_single_contig_1030181.start	225
    splits1 = nil
    splits2 = nil
    splits1 = row[0].match(/^(.+)\.(.+)$/)
    splits2 = row[1].match(/^(.+)\.(.+)$/)
    raise if splits1.nil? or splits2.nil?
    distance = row[2]
    row = [
      splits1[1], splits1[2], splits2[1], splits2[2], distance
      ].flatten
  end

  raise unless row.length == 5
  # e.g. seq1    end     seq23   start   6878

  probe1 = Probe.new
  probe1.contig_name = row[0]
  probe1.side = row[1]

  probe2 = Probe.new
  probe2.contig_name = row[2]
  probe2.side = row[3]

  if probe1.contig_name == probe2.contig_name and probe1.side == probe2.side
    num_probes_circular += 1
  else
    graph.add_edge probe1.to_setable, probe2.to_setable
  end
end

probes = graph.vertices.to_a

# Connect all starts to the ends
probes.each do |array|
  contig_name = array[0]

  start_probe = Probe.new
  start_probe.contig_name = contig_name
  start_probe.side = 'start'
  end_probe = Probe.new
  end_probe.contig_name = contig_name
  end_probe.side = 'end'

  graph.add_edge start_probe.to_setable, end_probe.to_setable
end
log.info "Removed #{num_probes_circular} connections that join a contig end to itself"

# First try the not computationally intensive way - can we find any?
edge_result = graph.some_edges_in_all_hamiltonian_cycles

cross_contig_connections = []
edge_result.edges_in_all.each do |v1, v2|
  if v1[0] != v2[0]
    cross_contig_connections.push [v1,v2]
  end
end

if !cross_contig_connections.empty?
  length = cross_contig_connections.length

  log.info "Good news. Found #{length} connections that are in all Hamiltonian paths (and thus can probably be scaffolded together):"
  cross_contig_connections.each do |connection|
    log.info connection[0].to_s + "\t" + connection[1].to_s
  end
  if length == probes.length and edge_result.contains_hamcycle != false
    log.info "Extra good news. You just scaffolded your genome into a single scaffold"
  end
end
if edge_result.contains_hamcycle == false
  log.warn "Bad news. The connectivity graph contains no Hamiltonian cycles, and so the contigs cannot be scaffolded into one circular genome"
end

# Determine if there are any ends that don't connect to anything
contig_names = []
Bio::FlatFile.foreach(options[:fasta_file]) do |seq|
  contig_name = seq.definition.split(/\s+/)[0]
  contig_names.push contig_name
  %w(start end).each do |side|
    probe = Probe.new
    probe.contig_name = contig_name
    probe.side = side
    if graph.edges[probe.to_setable].empty?
      log.info "Unable to find any connections from #{probe.to_setable}"
    end
  end
end

# Determine if there is any possible plasmids
num_plasmids_found = 0
contig_names.each do |contig_name|
  probe = Probe.new
  probe.contig_name = contig_name
  probe.side = 'start'
  rev_probe = Probe.new
  rev_probe.contig_name = contig_name
  rev_probe.side = 'end'
  if graph.edges[probe.to_setable].length == 1 and graph.edges[probe.to_setable].include?(rev_probe.to_setable)
    num_plasmids_found += 1
    log.info "Contig #{contig_name} appears to be circular and not connect to other contigs, suggesting it may be a plasmid"
  end
end
log.info "Found #{num_plasmids_found} contigs that appear to be plasmids based on connectivity"



log.info "Attempting a better but more computationally intensive method of determining edges that are in all hamiltonian paths.."
# First try to see if there is any hamiltonian paths?
paths = []
max_path_count = 4
operation_limit = 50000
graph.hamiltonian_cycles(operation_limit) do |path|
  if paths.length <= max_path_count
    paths.push path
  else
    break
  end
end

if paths.length < max_path_count
  log.info "Found exactly #{paths.length} Hamiltonian cycles"
else
  log.info "Gave up searching for Hamiltonian cycles as there are at least #{max_path_count} cycles"
end

# OK so
#edges_in_all = graph.edges_in_all_hamiltonian_cycles

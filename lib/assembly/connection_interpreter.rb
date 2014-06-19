require 'yargraph'


class Bio::FinishM::ConnectionInterpreter
  include Bio::FinishM::Logging

  # connections is an Enumerable of Probe object , sequences is a hash of name => DNA string
  def initialize(connections, sequences)
    @graph = Yargraph::UndirectedGraph.new
    @circular_probes = []
    @sequences = sequences
    @connections = connections

    # Add connections
    connections.each do |conn|
      if conn.probe1.to_settable == conn.probe2.to_settable
        @circular_probes.push con..probe1
      else
        @graph.add_edge conn.probe1.to_settable, conn.probe2.to_settable
      end
    end

    # Connect the start and the end of each probe together
    probes = graph.vertices.to_a
    sequences.each do |name, seq|
      start_probe = Probe.new
      start_probe.contig_name = name
      start_probe.side = :start
      end_probe = Probe.new
      end_probe.contig_name = name
      end_probe.side = :end

      graph.add_edge start_probe.to_settable, end_probe.to_settable
    end

    log.debug "Created a graph with #{@graph.vertices.to_a.length} verticies" if log.debug?
  end

  # Return sequences that exclusively connect the start to the end. In particular,
  # return the sequence name
  def circular_sequences
    to_return = []
    @connections.each do |conn|
      if conn.probe1.sequence_name == conn.probe2.sequence_name and
        conn.probe1.side != conn.probe2.side and
        @graph.edges[conn.probe1.to_settable].length == 1 and
        @graph.edges[conn.probe2.to_settable].length == 1

        to_return.push conn.probe1.sequence_name
      end
    end
    return to_return
  end


  # Return an Array of Connection objects that represent edges that must exist
  # in all hamiltonian cycles of the given contig set. These can likely
  # be scaffolded together.
  def likely_inter_contig_connections
    edge_result = @graph.some_edges_in_all_hamiltonian_cycles

    cross_contig_connections = []
    edge_result.edges_in_all.each do |v1, v2|
      if v1[0] != v2[0] #dont count edges that connect starts to ends of the same contig
        conn = Connection.new
        conn.probe1 = Probe.new(v1)
        conn.probe2 = Probe.new(v2)
        cross_contig_connections.push conn
      end
    end

    return cross_contig_connections
  end

  # Single linkage cluster the likely_inter_contig_connections
  # and the start to ends for each of the contigs.
  def scaffolds
    # It is like an
    # assembly problem because each vertex can only be connected to
    # two others - 1 intra-contig and 1 inter-contig (unless it is circular)
    likelies = likely_inter_contig_connections
    likelies_edge_set = EdgeSet.new
    likelies.each do |conn|
      likelies_edge_set.add_edge conn.probe1.to_settable, conn.probe2.to_settable
    end

    scaffolded_paths = path

    # while there is more elements in the likelies set
    while !likelies.empty?
      # 'pop' an arbitrary edge out
      starting_edge = likelies_edge_set.to_a[0]
      likelies_edge_set.delete(starting_edge[0], starting_edge[1])

      # go 'left'. Connect the other side of the left.
      lefts = [starting_edge[0]]
      # while there is another node to the left
      while next_edge = likelies_edge_set[lefts[-1]]
        likelies_edge_set.delete next_edge[0], next_edge[1]
        lefts.push next_edge
      end
      # and go right
      rights = [starting_edge[1]]
      while next_edge = likelies_edge_set[rights[-1]]
        likelies_edge_set.delete next_edge[0], next_edge[1]
        lefts.push next_edge
      end

      # Add the left and the right together into one path
      path = left.reverse+right[1..-1]
    end

    # for each scaffolded set, create new scaffold object
    scaffolds = []
    scaffolded_contigs = Set.new
    scaffolded_paths.each do |path|
      raise if path.length % 2 != 0
      scaffold = Bio::FinishM::ScaffoldBreaker::Scaffold.new
      scaff.name = "scaffold#{scaffolds.length+1}"
      path.each_with_index do |probe, i|
        next if i % 2 == 1
        contig = Bio::FinishM::ScaffoldBreaker::UnscaffoldedContig.new
        contig.scaffold = scaffold
        seq = @sequences[probe.sequence_name]
        if probe.first_side == :start
          contig.sequence = seq
        else
          contig.sequence = Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
        end
        contig.scaffold_position_start = (scaffold.sequence+('N'*10)).length
        contig.scaffold_position_end = contig.scaffold_position_start + seq.length - 1
        scaffold.contigs.push contig
      end
      scaffolds.push scaffold
    end

    # for each contig that is not in a contig, add as singleton
    @sequences.each do |name, seq|
      unless scaffolded_contigs.include?(name)
        scaff = Bio::FinishM::ScaffoldBreaker::Scaffold.new
        scaff.name = "scaffold#{scaffolds.length+1}"
        contig = Bio::FinishM::ScaffoldBreaker::UnscaffoldedContig.new
        contig.scaffold = scaff
        contig.sequence = seq
        contig.scaffold_position_start = 1
        contig.scaffold_position_end = seq.length
        scaffolds.push scaff
      end
    end

    return scaffolds
  end

  class Connection
    attr_accessor :probe1, :probe2

    attr_accessor :distance
  end

  class Probe
    attr_accessor :side #:start or :end
    attr_accessor :sequence_name #name of the underlying sequence

    def initialize(settable_representation=nil)
      unless settable_representation.nil?
        @sequence_name = settable_representation[0]
        @side = settable_representation[1]
      end
    end

    def to_settable
      [@sequence_name, @side]
    end
  end
end
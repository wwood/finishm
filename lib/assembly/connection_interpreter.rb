require '/home/ben/git/yargraph/lib/yargraph'


class Bio::FinishM::ConnectionInterpreter
  include Bio::FinishM::Logging

  # connections is an Enumerable of Probe object , sequences is a hash of name => DNA string
  def initialize(connections, sequences)
    @graph = Yargraph::UndirectedGraph.new
    @circular_probes = []
    @sequences = sequences

    # Setup hash of setable to original
    # Assume there is only 1 connection between two contig ends
    @connection_hash = {}
    connections.each do |conn|
      key = conn.to_settable
      raise "Duplicate connections not handled (yet?), found #{conn} => #{key}" if @connection_hash.key?(key)
      @connection_hash[key] = conn
    end

    # Add connections
    connections.each do |conn|
      if conn.probe1.to_settable == conn.probe2.to_settable
        @circular_probes.push con..probe1
      else
        @graph.add_edge conn.probe1.to_settable, conn.probe2.to_settable
      end
    end

    log.debug "Created a graph with #{@graph.vertices.to_a.length} vertices and #{@graph.edges.length} edges" if log.debug?
  end

  def connections
    @connection_hash.values
  end

  # Return sequences that exclusively connect the start to the end. In particular,
  # return an Array of sequence names
  def circular_sequences
    to_return = []
    connections.each do |conn|
      if conn.probe1.sequence_name == conn.probe2.sequence_name and
        conn.probe1.side != conn.probe2.side and
        @graph.edges[conn.probe1.to_settable].length == 1 and
        @graph.edges[conn.probe2.to_settable].length == 1

        to_return.push conn.probe1.sequence_name
      end
    end
    return to_return
  end


  # Return an Array of Connection objects that represent edges where
  # there is only a single connection from both side
  def doubly_single_contig_connections
    likelies = []

    already_seen_connections = Set.new

    @graph.vertices.each do |v|
      # If there is only 1 connection on both sides, then go with that
      neighbours = @graph.neighbours(v)
      log.debug "Testing connection between #{v} and #{neighbours}"
      if neighbours.length == 1 and @graph.neighbours(neighbours[0]).length == 1
        log.debug "Connection passed the doubly-test" if log.debug?
        neighbour = neighbours[0]

        conn = Connection.new
        conn.probe1 = Probe.new(v)
        conn.probe2 = Probe.new(neighbour)
        settable = conn.to_settable
        # Record the connection unless it is duplicate
        unless already_seen_connections.include?(settable)
          likelies.push @connection_hash[settable]
          already_seen_connections << settable
        end
      end
    end

    return likelies
  end

  # Single linkage cluster the likely_inter_contig_connections
  # and the start to ends for each of the contigs. Assumes
  def scaffolds(contig_connections)
    # It is like an (easy)
    # assembly problem because each vertex can only be connected to
    # two others - 1 intra-contig and 1 inter-contig (unless it is circular)
    contig_connections ||= doubly_single_contig_connections
    likelies_edge_set = Yargraph::UndirectedGraph::EdgeSet.new
    contig_connections.each do |conn|
      likelies_edge_set.add_edge conn.probe1.to_settable, conn.probe2.to_settable
    end

    scaffolded_paths = []
    circular_single_contigs = Set.new

    # while there is more elements in the likelies set,
    # 'pop' an arbitrary edge out of the graph
    while starting_edge = likelies_edge_set.pop
      log.debug "starting to scaffold from #{starting_edge}" if log.debug?

      # Ignore likelies that are circular
      if starting_edge[0][0] == starting_edge[1][0]
        log.debug "Not scaffolding contig #{starting_edge[0][0] } since it appears to be circular" if log.debug?
        circular_single_contigs << starting_edge[0][0]
        next
      end

      circular = false

      # go 'left'. Connect the other side of the left.
      lefts = [Probe.new(starting_edge[0])]
      rights = [Probe.new(starting_edge[1])]
      log.debug "rights was #{rights[0].to_s}" if log.debug?
      # while there is another node to the left
      while next_probe = likelies_edge_set[lefts[-1].companion.to_settable].to_a[0]
        next_probe_probe = Probe.new(next_probe)
        companion = lefts[-1].companion

        likelies_edge_set.delete next_probe, companion.to_settable
        if next_probe_probe.companion.to_settable == rights[0].to_settable
          log.debug "Found multi-contig circularity between #{next_probe_probe.companion} and #{rights[0]}" if log.debug?
          circular = true
          break
        end

        lefts.push companion
        lefts.push next_probe_probe
        log.debug "Adding node to the left: #{next_probe} and companion #{companion}" if log.debug?
      end
      # and go right
      while next_probe = likelies_edge_set[rights[-1].companion.to_settable].to_a[0]
        companion = rights[-1].companion
        rights.push companion
        rights.push Probe.new(next_probe)
        log.debug "Adding node to the right: #{next_probe} and companion #{companion}" if log.debug?
        likelies_edge_set.delete next_probe, companion.to_settable
      end

      # Add the left and the right together into one path
      scaffolded_paths.push(
        PossiblyCircularArray.new(
          [lefts[-1].companion]+
            lefts.reverse+
            rights+
            [rights[-1].companion],
          circular)
        )
    end
    if log.debug?
      log.debug "Found #{scaffolded_paths.length} multi-contig scaffold(s):"
      scaffolded_paths.each do |path|
        log.debug "Scaffold: #{path.collect{|e| e.to_s}.join(', ')}"
      end
    end


    # for each scaffolded set, create new scaffold object
    scaffolds = []
    scaffolded_contigs = Set.new
    scaffolded_paths.each do |path|
      raise if path.length % 2 != 0
      scaffold = Scaffold.new
      scaffold.circular = path.circular
      scaffold.name = "scaffold#{scaffolds.length+1}"
      path.each_with_index do |probe, i|
        next if i % 2 == 1
        contig = UnscaffoldedContig.new
        contig.original_name = probe.sequence_name
        contig.scaffold = scaffold
        seq = @sequences[probe.sequence_name]
        if probe.side == :start
          contig.sequence = seq
        else
          contig.sequence = Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
        end
        scaffold.contigs ||= []
        contig.scaffold_position_start = scaffold.sequence.length
        contig.scaffold_position_start += 10 unless scaffold.contigs.empty?
        contig.scaffold_position_end = contig.scaffold_position_start + seq.length - 1
        scaffold.contigs.push contig
        scaffolded_contigs << probe.sequence_name
      end
      scaffolds.push scaffold
    end

    # for each contig that is not in a contig, add as singleton
    @sequences.each do |name, seq|
      unless scaffolded_contigs.include?(name)
        scaff = Scaffold.new
        scaff.name = "singleton#{scaffolds.length+1}"
        contig = UnscaffoldedContig.new
        contig.original_name = name
        contig.scaffold = scaff
        contig.sequence = seq
        contig.scaffold_position_start = 0
        contig.scaffold_position_end = seq.length-1
        scaff.contigs = [contig]
        if circular_single_contigs.include?(name)
          scaff.circular = true
        else
          scaff.circular = false
        end
        scaffolds.push scaff
      end
    end

    return scaffolds
  end

  class Connection
    # Probe objects
    attr_accessor :probe1, :probe2

    attr_accessor :distance

    def to_s
      [@probe1, @probe2].join('/')+":#{@distance}"
    end

    def to_settable
      if @probe1.sequence_name < @probe2.sequence_name
        return [@probe1.to_settable, @probe2.to_settable].flatten
      elsif @probe1.sequence_name == @probe2.sequence_name
        if @probe1.side < @probe2.side
          return [@probe1.to_settable, @probe2.to_settable].flatten
        else
          return [@probe2.to_settable, @probe1.to_settable].flatten
        end
      else
        return [@probe2.to_settable, @probe1.to_settable].flatten
      end
    end
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

    def to_s
      side = @side == :start ? 's' : 'e'
      "#{@sequence_name}#{side}"
    end

    # Return a probe representing the other side of the contig
    def companion
      companion = Probe.new
      companion.sequence_name = @sequence_name
      companion.side = @side == :start ? :end : :start
      return companion
    end
  end

  class Scaffold < Bio::FinishM::ScaffoldBreaker::Scaffold
    attr_accessor :circular
    def circular?
      @circular
    end
  end

  class UnscaffoldedContig < Bio::FinishM::ScaffoldBreaker::UnscaffoldedContig
    attr_accessor :original_name
  end

  class PossiblyCircularArray < Array
    attr_accessor :circular

    def initialize(array, circular)
      @circular = circular
      super(array)
    end
  end
end

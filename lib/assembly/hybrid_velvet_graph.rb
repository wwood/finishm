# A velvet Graph where the nodes and arcs are in Ruby, but the NodedReads are in C
class Bio::FinishM::HybridGraph < Bio::Velvet::Graph
  include Bio::FinishM::Logging

  attr_accessor :bio_velvet_underground_graph

  def initialize(bio_velvet_graph, bio_velvet_underground_graph)
    @bio_velvet_underground_graph = bio_velvet_underground_graph

    @nodes = NodeArray.new(bio_velvet_graph, bio_velvet_underground_graph, self)
    @hash_length = bio_velvet_graph.hash_length
    @number_of_nodes = bio_velvet_graph.number_of_nodes
    @number_of_sequences = bio_velvet_graph.number_of_sequences
    @arcs = bio_velvet_graph.arcs
  end

  class NodeArray
    include Enumerable

    def initialize(bio_velvet_graph, bio_velvet_underground_graph, parent_graph)
      @bio_velvet_graph = bio_velvet_graph
      @bio_velvet_underground_graph = bio_velvet_underground_graph
      @parent_graph = parent_graph
    end

    def []=(node_id, value)
      raise "method not implemented"
    end

    def [](node_id)
      bio_velvet_node = @bio_velvet_graph.nodes[node_id]
      bio_velvet_node.short_reads = LazyShortReadArray.new(@bio_velvet_underground_graph.nodes[node_id])
      bio_velvet_node.parent_graph = @parent_graph
      return bio_velvet_node
    end

    def delete(node)
      @bio_velvet_graph.nodes.delete node
    end

    def length
      @bio_velvet_graph.nodes.length
    end

    def each(&block)
      @bio_velvet_graph.nodes.each do |node|
        block.yield self[node.node_id]
      end
    end
  end

  # Hold pointers to short reads, but don't actually create them until required.
  # Utilizes method_missing to pick up calls made to the underlying object
  class LazyShortReadArray
    @short_reads = nil

    def method_missing(m, *args, &block)
      read_short_reads.send(m, *args, &block)
    end

    def read_short_reads
      @short_reads ||= @parent_node.short_reads
    end

    def initialize(node)
      @parent_node = node
    end
  end
end


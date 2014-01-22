$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rspec'
require 'priner'

# Requires supporting files with custom matchers and macros, etc,
# in ./support/ and its subdirectories.
Dir["#{File.dirname(__FILE__)}/support/**/*.rb"].each {|f| require f}

RSpec.configure do |config|

end
TEST_DATA_DIR = File.join(File.dirname(__FILE__),'data')
FINISHM_SCRIPT_PATH = File.join(File.dirname(__FILE__),'..','bin','finishm')

class GraphTesting
  def self.emit(arc_pairs)
    node_id_to_node = {}
    graph = Bio::Velvet::Graph.new
    arc_array = Bio::Velvet::Graph::ArcArray.new
    graph.arcs = arc_array
    nodes = Bio::Velvet::Graph::NodeArray.new
    graph.nodes = nodes
    graph.hash_length = 7

    arc_pairs.each do |node_ids|
      raise unless node_ids.length == 2

      # Create the nodes if necessary
      node_ids.each_with_index do |ident|
        node = node_id_to_node[ident]
        if node.nil?
          node = Bio::Velvet::Graph::Node.new
          node.node_id = ident
          node_id_to_node[ident] = node
          node.parent_graph = graph
          node.ends_of_kmers_of_node = 'A'*10
          node.ends_of_kmers_of_twin_node = 'T'*10
        end
        nodes[ident] = node
      end

      # Create the arc itself
      arc = Bio::Velvet::Graph::Arc.new
      arc.begin_node_id = node_ids[0]
      arc.end_node_id = node_ids[1]
      arc.begin_node_direction = true
      arc.end_node_direction = true
      arc_array.push arc
    end

    return graph
  end

  def self.emit_ss(arc_pairs, start_node_id, stop_node_id)
    graph = emit(arc_pairs)
    initial = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    initial.node = graph.nodes[start_node_id]
    initial.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
    initial_path.add_oriented_node initial
    terminal = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    terminal.node = graph.nodes[stop_node_id]
    terminal.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

    return graph, initial_path, terminal
  end

  def self.finishm_graph(arc_pairs, probes)
    graph = emit(arc_pairs)
    finishm_graph = Bio::FinishM::ProbedGraph.new
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = probes.collect{|probe| graph.nodes[probe]}
    finishm_graph.probe_node_directions = probes.collect{true}

    return finishm_graph
  end

  def self.sorted_paths(paths)
    paths.collect do |path|
      path.collect{|n| n.node_id}
    end.sort
  end

  def self.sorted_array_of_paths(path_array)
    path_array.collect do |paths|
      paths.collect do |path|
        path.collect{|n| n.node_id}
      end.sort
    end
  end

  def self.sorted_fluffers(path_array)
    path_array.collect do |fluffer_path_set|
      paths_and_fates = []
      fluffer_path_set.each_with_index do |path, i|
        paths_and_fates.push [
          path.collect{|n| n.node_id},
          fluffer_path_set.fates[i]
        ]
      end
      paths_and_fates.sort
    end
  end
end

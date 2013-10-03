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

class GraphTesting
  def self.emit(arc_pairs)
    node_id_to_node = {}
    graph = Bio::Velvet::Graph.new
    arc_array = Bio::Velvet::Graph::ArcArray.new
    graph.arcs = arc_array
    nodes = Bio::Velvet::Graph::NodeArray.new
    graph.nodes = nodes

    arc_pairs.each do |node_ids|
      raise unless node_ids.length == 2

      # Create the nodes if necessary
      node_ids.each_with_index do |ident|
        node = node_id_to_node[ident]
        if node.nil?
          node = Bio::Velvet::Graph::Node.new
          node.node_id = ident
          node_id_to_node[ident] = node
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

  def self.sorted_paths(paths)
    paths.collect do |path|
      path.collect{|n| n.node_id}
    end.sort
  end
end

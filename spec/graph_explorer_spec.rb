require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
class Util
  def self.revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "GraphExplorer" do

  it 'should find one easy path' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
    ])
    initial_node = graph.nodes[1]
    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    initial_trail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    paths = Bio::AssemblyGraphAlgorithms::GraphExplorer.new.explore_from_node(graph, initial_trail, nil)
    paths.collect{|path| path.termination_type}.should == [
      'Dead end / coverage',
    ]
    GraphTesting.sorted_paths(paths.collect{|path| path.path}).should == [
      [1,2,3,4],
    ]
  end

  it 'should find two easy paths' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [3,5],
    ])
    initial_node = graph.nodes[1]
    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    initial_trail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    paths = Bio::AssemblyGraphAlgorithms::GraphExplorer.new.explore_from_node(graph, initial_trail, nil)
    paths.collect{|path| path.termination_type}.should == [
      'Dead end / coverage',
      'Dead end / coverage',
    ]
    GraphTesting.sorted_paths(paths.collect{|path| path.path}).should == [
      [1,2,3,4],
      [1,2,3,5],
    ]
  end

  it 'should find two paths with bubble' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [3,5],
      [4,6],
      [5,6],
      [6,7],
    ])
    initial_node = graph.nodes[1]
    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    initial_trail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    paths = Bio::AssemblyGraphAlgorithms::GraphExplorer.new.explore_from_node(graph, initial_trail, nil)
    paths.collect{|path| path.termination_type}.should == [
      'Dead end / coverage',
      'Dead end / coverage',
    ]
    GraphTesting.sorted_paths(paths.collect{|path| path.path}).should == [
      [1,2,3,4,6,7],
      [1,2,3,5,6,7],
    ]
  end

  it 'should handle loops' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,1],
      [3,4],
    ])
    initial_node = graph.nodes[1]
    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    initial_trail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    paths = Bio::AssemblyGraphAlgorithms::GraphExplorer.new.explore_from_node(graph, initial_trail, nil)
    GraphTesting.sorted_paths(paths.collect{|path| path.path}).should == [
      [1,2,3,1],
      [1,2,3,4],
    ]
    paths.collect{|path| path.termination_type}.should == [
      'Loop',
      'Dead end / coverage',
    ]
  end

#  it 'should find paths not both ending at terminal node' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [1,5],
#      [5,3]
#    ])
#    initial_node = graph.nodes[1]
#    initial_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
#    initial_trail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail.START_IS_FIRST
#    cartographer = Bio::AssemblyGraphAlgorithms::GraphExplorer.new
#    paths = cartographer.explore_from_node(graph, initial_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,4],
#      [1,5,3,4],
#    ]
#  end
#
#  it 'should find through consecutive loops ending at terminal' do
#    # 1 2/3 4 5/6 7
#    graph = GraphTesting.emit([
#      [1,2],
#      [1,3],
#      [2,4],
#      [3,4],
#      [4,5],
#      [4,6],
#      [5,7],
#      [6,7],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[7]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,4,5,7],
#      [1,2,4,6,7],
#      [1,3,4,5,7],
#      [1,3,4,6,7],
#    ]
#  end
#
#    it 'should find through consecutive loops not ending at terminal' do
#    # 1 2/3 4 5/6 7 8
#    graph = GraphTesting.emit([
#      [1,2],
#      [1,3],
#      [2,4],
#      [3,4],
#      [4,5],
#      [4,6],
#      [5,7],
#      [6,7],
#      [7,8]
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[8]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,4,5,7,8],
#      [1,2,4,6,7,8],
#      [1,3,4,5,7,8],
#      [1,3,4,6,7,8],
#    ]
#  end
#
#  it 'should find loop off loop' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [2,4],
#      [3,8],
#      [4,5],
#      [4,6],
#      [5,7],
#      [6,7],
#      [7,8],
#      [8,9],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[9]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,8,9],
#      [1,2,4,5,7,8,9],
#      [1,2,4,6,7,8,9],
#    ].sort
#  end
#
#  it 'should find loop off loop with three way' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [1,3],
#      [1,4],
#      [2,5],
#      [3,5],
#      [4,8],
#      [5,6],
#      [5,7],
#      [6,9],
#      [7,9],
#      [8,9],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[9]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,5,6,9],
#      [1,2,5,7,9],
#      [1,3,5,6,9],
#      [1,3,5,7,9],
#      [1,4,8,9],
#    ].sort
#  end
#
#  it 'should not fail when there is one path beyond the leash (1st), and another not (2nd)' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [1,3],
#      [3,4],
#    ])
#    graph.hash_length = 87
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[4]
#    graph.nodes[1].ends_of_kmers_of_node = 'A'*10
#    graph.nodes[2].ends_of_kmers_of_node = 'A'*100
#    graph.nodes[3].ends_of_kmers_of_node = 'A'*10
#    graph.nodes[4].ends_of_kmers_of_node = 'A'*10
#    (1..4).each do |node_id|
#      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#    end
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 200, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,4],
#      [1,3,4],
#    ].sort
#
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 100, true)
#    (1..4).each do |node_id|
#      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#    end
#    GraphTesting.sorted_paths(paths).should == [
#      [1,3,4],
#    ].sort
#
#  end
#
#
#  it 'probably fails with the nasty leash bug, which is hard to fix' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [4,5],
#      [1,3],
#    ])
#    graph.hash_length = 87
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[5]
#    graph.nodes[1].ends_of_kmers_of_node = 'A'*10
#    graph.nodes[2].ends_of_kmers_of_node = 'A'*70
#    graph.nodes[3].ends_of_kmers_of_node = 'A'*15
#    graph.nodes[4].ends_of_kmers_of_node = 'A'*10
#    graph.nodes[5].ends_of_kmers_of_node = 'A'*10
#    (1..4).each do |node_id|
#      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#    end
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 200, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,4,5],
#      [1,3,4,5],
#    ].sort
#
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 100, true)
#    (1..4).each do |node_id|
#      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#    end
#    GraphTesting.sorted_paths(paths).should == [
#      [1,3,4,5],
#    ].sort
#
#  end
#
#  it 'should not fail in this special case I realised might trip up the algorithm' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [4,5],
#
#      [1,6],
#      [6,3],
#      [6,5],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[5]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,4,5],
#      [1,2,3,4,6,5],
#      [1,6,3,4,5],
#      [1,6,5],
#    ].sort
#  end
#
#  it 'should not fail in another special case I realised might trip up the algorithm' do
#    #NOTE: to fix this, one must first fix the above graph problem. Argh.
#    # The problem is that a simple joining of the golden path 1,2,3,4,5 and the
#    # golden fragment 1,6,3 yields a circular path
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [4,5],
#
#      [1,6],
#      [6,3],
#      [6,5],
#
#      [4,6],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[5]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3,4,5],
#      [1,2,3,4,6,5],
#      [1,6,3,4,5],
#      [1,6,5],
#    ].sort
#  end
#
#  it 'should give the same answer as a more straightfoward (naive?) repeated depth first search style' do
#    raise "need to do some simulation work here to write the test"
#  end
#
#  it 'should not get confused by a 1 node cycle' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,2],
#      [2,3],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[3]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3],
#    ].sort
#  end
#
#  it 'should not get confused by a 2 node cycle' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,4],
#      [4,2],
#      [2,3],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[3]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
#    GraphTesting.sorted_paths(paths).should == [
#      [1,2,3],
#    ].sort
#  end
#
#  it 'should give the correct part1 for a 2 node cycle' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,4],
#      [4,2],
#      [2,3],
#    ])
#    initial_node = graph.nodes[1]
#    terminal_node = graph.nodes[3]
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#
#    initial_path = Bio::Velvet::Graph::OrientedNodeTrail.new
#    way = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
#    initial_path.add_node initial_node, way
#
#    half_result = cartographer.find_all_trails_squid_way_part1(graph, initial_path, terminal_node, nil)
#    half_result.golden_path.collect{|n| n.node.node_id}.should == [1,2,3]
#    half_result.golden_fragments.collect{|fragment| fragment.collect{|n| n.node.node_id}}.should == []
#  end
#
#  it 'should calculate connections using a simple depth first search' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#    ])
#    finishm_graph = Bio::FinishM::ProbedGraph.new
#    finishm_graph.graph = graph
#    finishm_graph.probe_nodes = [graph.nodes[1],graph.nodes[2]]
#    finishm_graph.probe_node_directions = [true, false]
#
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
#      [0,1] => 10
#    }
#  end
#
#  it 'should calculate connections using a simple depth first search multiple singly connected' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [4,5]
#    ])
#    finishm_graph = Bio::FinishM::ProbedGraph.new
#    finishm_graph.graph = graph
#    finishm_graph.probe_nodes = [
#      graph.nodes[1],graph.nodes[3],
#      graph.nodes[4],graph.nodes[5],
#    ]
#    finishm_graph.probe_node_directions = [true, false, true, false]
#
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
#      [0,1] => 20,
#      [2,3] => 10,
#    }
#  end
#
#  it 'should calculate connections using a simple depth first search multiply connected' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [4,5],
#    ])
#    finishm_graph = Bio::FinishM::ProbedGraph.new
#    finishm_graph.graph = graph
#    finishm_graph.probe_nodes = [
#      graph.nodes[1],graph.nodes[3],
#      graph.nodes[4],graph.nodes[5],
#    ]
#    finishm_graph.probe_node_directions = [true, false, true, false]
#
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
#      [0,1] => 20,
#      [0,3] => 40,
#      [2,3] => 10,
#    }
#  end
#
#  it 'should calculate connections using a simple depth first search respect leash' do
#    graph = GraphTesting.emit([
#      [1,2],
#      [2,3],
#      [3,4],
#      [4,5],
#    ])
#    finishm_graph = Bio::FinishM::ProbedGraph.new
#    finishm_graph.graph = graph
#    finishm_graph.probe_nodes = [
#      graph.nodes[1],graph.nodes[3],
#      graph.nodes[4],graph.nodes[5],
#    ]
#    finishm_graph.probe_node_directions = [true, false, true, false]
#
#    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
#    cartographer.depth_first_search_with_leash(finishm_graph, 20).should == {
#      [0,1] => 20,
#      #[0,3] => 40,
#      [2,3] => 10,
#    }
#  end
end

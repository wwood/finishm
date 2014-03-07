require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
class Util
  def self.revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::LoggerPlus.new('bio-velvet'); Bio::Log::CLI.configure('bio-velvet')

describe "AcyclicConnectionFinder" do

  it 'should calculate an easy trail' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    # Make sure we are on the same page
    graph.arcs.length.should eq(4)
    graph.nodes.length.should eq(4)
    nodes = graph.nodes

    #viser = Bio::Assembly::ABVisualiser.new
    #gv = viser.graphviz(graph)
    #gv.output :svg => 'a.svg', :use => :neato


    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new

    log = nil
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [1],
    ]

    log.debug "============= starting second test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[2], true)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    log.debug "Found trails: #{trails[0].collect{|node| node.node_id}.join(',')}" if log
    GraphTesting.sorted_paths(trails).should == [
      [1,2],
    ]

    log.debug "============= starting third test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], true)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [1,3],
    ]

    log.debug "============= starting 4th test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], false)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [1,2,4],
      [1,3,4],
    ].sort
  end

  it 'should calculate trails in reverse' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails_reverse','Assem','LastGraph')
    # Make sure we are on the same page
    graph.arcs.length.should eq(4)
    graph.nodes.length.should eq(4)
    nodes = graph.nodes

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new

    log = nil
    #Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [1],
    ]

    log.debug "============= starting second test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], false)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[2], false)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    log.debug "Found trails: #{trails[0].collect{|node| node.node_id}.join(',')}" if log
    GraphTesting.sorted_paths(trails).should == [
      [3,2],
    ]

    log.debug "============= starting third test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], false)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], false)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [3,4],
    ]

    log.debug "============= starting 4th test" if log
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], false)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], false)
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [3,2,1],
      [3,4,1],
    ]
  end

  it 'should not go overboard searching for paths' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    # Make sure we are on the same page
    graph.arcs.length.should eq(4)
    graph.nodes.length.should eq(4)
    nodes = graph.nodes

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new

    log = nil
    #Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], false)

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [[3,4]]

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 20)
    GraphTesting.sorted_paths(trails).should == []


    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], false)

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999)
    GraphTesting.sorted_paths(trails).should == [
      [1,2,4],
      [1,3,4],
    ]

    # Ignoring the below test because the squid way doesn't violate the leash to get both paths,
    # where the old method did.
    #    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 228+29+224+30+2, true)
    #    GraphTesting.sorted_paths(trails).should == [
    #      [1,2,4]
    #    ]
  end

  it 'should find paths not both ending at terminal node' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [1,5],
      [5,3]
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,4],
      [1,5,3,4],
    ]
  end

  it 'should find through consecutive loops ending at terminal' do
    # 1 2/3 4 5/6 7
    graph = GraphTesting.emit([
      [1,2],
      [1,3],
      [2,4],
      [3,4],
      [4,5],
      [4,6],
      [5,7],
      [6,7],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[7], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,4,5,7],
      [1,2,4,6,7],
      [1,3,4,5,7],
      [1,3,4,6,7],
    ]
  end

    it 'should find through consecutive loops not ending at terminal' do
    # 1 2/3 4 5/6 7 8
    graph = GraphTesting.emit([
      [1,2],
      [1,3],
      [2,4],
      [3,4],
      [4,5],
      [4,6],
      [5,7],
      [6,7],
      [7,8]
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[8], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,4,5,7,8],
      [1,2,4,6,7,8],
      [1,3,4,5,7,8],
      [1,3,4,6,7,8],
    ]
  end

  it 'should find loop off loop' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [2,4],
      [3,8],
      [4,5],
      [4,6],
      [5,7],
      [6,7],
      [7,8],
      [8,9],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[9], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,8,9],
      [1,2,4,5,7,8,9],
      [1,2,4,6,7,8,9],
    ].sort
  end

  it 'should find loop off loop with three way' do
    graph = GraphTesting.emit([
      [1,2],
      [1,3],
      [1,4],
      [2,5],
      [3,5],
      [4,8],
      [5,6],
      [5,7],
      [6,9],
      [7,9],
      [8,9],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[9], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,5,6,9],
      [1,2,5,7,9],
      [1,3,5,6,9],
      [1,3,5,7,9],
      [1,4,8,9],
    ].sort
  end

  it 'should not fail when there is one path beyond the leash (1st), and another not (2nd)' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [1,3],
      [3,4],
    ])
    graph.hash_length = 87
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[4], true)
    graph.nodes[1].ends_of_kmers_of_node = 'A'*10
    graph.nodes[2].ends_of_kmers_of_node = 'A'*100
    graph.nodes[3].ends_of_kmers_of_node = 'A'*10
    graph.nodes[4].ends_of_kmers_of_node = 'A'*10
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 200)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,4],
      [1,3,4],
    ].sort

    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 100)
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    GraphTesting.sorted_paths(paths).should == [
      [1,3,4],
    ].sort

  end


  it 'probably fails with the nasty leash bug, which is hard to fix' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [4,5],
      [1,3],
    ])
    graph.hash_length = 87
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[5], true)
    graph.nodes[1].ends_of_kmers_of_node = 'A'*10
    graph.nodes[2].ends_of_kmers_of_node = 'A'*70
    graph.nodes[3].ends_of_kmers_of_node = 'A'*15
    graph.nodes[4].ends_of_kmers_of_node = 'A'*10
    graph.nodes[5].ends_of_kmers_of_node = 'A'*10
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 200)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,4,5],
      [1,3,4,5],
    ].sort

    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 60)
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    GraphTesting.sorted_paths(paths).should == [
      [1,3,4,5],
    ].sort

  end

  it 'should not fail in this special case I realised might trip up the algorithm' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [4,5],

      [1,6],
      [6,3],
      [6,5],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[5], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,4,5],
      [1,6,3,4,5],
      [1,6,5],
    ].sort
  end

  it 'should give the same answer as a more straightfoward (naive?) repeated depth first search style' do
    raise "need to do some simulation work here to write the test"
  end

  it 'should not get confused by a 1 node cycle' do
    graph = GraphTesting.emit([
      [1,2],
      [2,2],
      [2,3],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths).should == [
    ].sort
  end

  it 'should not get confused by a 2 node cycle' do
    graph = GraphTesting.emit([
      [1,2],
      [2,4],
      [4,2],
      [2,3],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[3], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil)
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths).should == [
    ].sort
  end

  it 'should calculate connections using a simple depth first search' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
    ])
    finishm_graph = Bio::FinishM::ProbedGraph.new
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = [graph.nodes[1],graph.nodes[2]]
    finishm_graph.probe_node_directions = [true, false]

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
      [0,1] => 10
    }
  end

  it 'should calculate connections using a simple depth first search multiple singly connected' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [4,5]
    ])
    finishm_graph = Bio::FinishM::ProbedGraph.new
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = [
      graph.nodes[1],graph.nodes[3],
      graph.nodes[4],graph.nodes[5],
    ]
    finishm_graph.probe_node_directions = [true, false, true, false]

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
      [0,1] => 20,
      [2,3] => 10,
    }
  end

  it 'should calculate connections using a simple depth first search multiply connected' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [4,5],
    ])
    finishm_graph = Bio::FinishM::ProbedGraph.new
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = [
      graph.nodes[1],graph.nodes[3],
      graph.nodes[4],graph.nodes[5],
    ]
    finishm_graph.probe_node_directions = [true, false, true, false]

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    cartographer.depth_first_search_with_leash(finishm_graph, 10000).should == {
      [0,1] => 20,
      [0,3] => 40,
      [2,3] => 10,
    }
  end

  it 'should calculate connections using a simple depth first search respect leash' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4],
      [4,5],
    ])
    finishm_graph = Bio::FinishM::ProbedGraph.new
    finishm_graph.graph = graph
    finishm_graph.probe_nodes = [
      graph.nodes[1],graph.nodes[3],
      graph.nodes[4],graph.nodes[5],
    ]
    finishm_graph.probe_node_directions = [true, false, true, false]

    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    cartographer.depth_first_search_with_leash(finishm_graph, 20).should == {
      [0,1] => 20,
      #[0,3] => 40,
      [2,3] => 10,
    }
  end

it 'should find a path when the initial is the terminal and directions match' do
    graph = GraphTesting.emit([
      [1,2],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999)
    GraphTesting.sorted_paths(paths).should == [
      [1],
      ].sort
  end

  it 'should not find a path when the initial is the terminal but the directions do not match' do
    graph = GraphTesting.emit([
      [1,2],
    ])
    initial_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
    terminal_node = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], false)
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999)
    GraphTesting.sorted_paths(paths).should == [
      ].sort
  end
end

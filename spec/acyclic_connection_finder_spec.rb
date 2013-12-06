require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
class Util
  def self.revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[1]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 1
    trails[0].length.should == 1
    trails[0][0].node_id.should == 1

    log.debug "============= starting second test" if log
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[2]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 1
    log.debug "Found trails: #{trails[0].collect{|node| node.node_id}.join(',')}" if log
    trails[0].length.should == 2
    trails[0][0].node_id.should == 1
    trails[0][1].node_id.should == 2

    log.debug "============= starting third test" if log
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[3]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 1
    trails[0].length.should == 2
    trails[0][0].node_id.should == 1
    trails[0][1].node_id.should == 3

    log.debug "============= starting 4th test" if log
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[4]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 2
    trails[0].length.should == 3
    trails[0][0].node_id.should == 1
    [2,3].include?(trails[0][1].node_id).should == true
    trails[0][2].node_id.should == 4
    trails[1][0].node_id.should == 1
    [2,3].include?(trails[1][1].node_id).should == true
    trails[1][1].node_id.should_not == trails[0][1].node_id
    trails[1][2].node_id.should == 4
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[1]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 1
    trails[0].length.should == 1
    trails[0][0].node_id.should == 1

    log.debug "============= starting second test" if log
    initial_node = graph.nodes[3]
    terminal_node = graph.nodes[2]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, false)
    trails.length.should == 1
    log.debug "Found trails: #{trails[0].collect{|node| node.node_id}.join(',')}" if log
    trails[0].length.should == 2
    trails[0][0].node_id.should == 3
    trails[0][1].node_id.should == 2

    log.debug "============= starting third test" if log
    initial_node = graph.nodes[3]
    terminal_node = graph.nodes[4]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, false)
    trails.length.should == 1
    trails[0].length.should == 2
    trails[0][0].node_id.should == 3
    trails[0][1].node_id.should == 4

    log.debug "============= starting 4th test" if log
    initial_node = graph.nodes[3]
    terminal_node = graph.nodes[1]
    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, false)
    trails.length.should == 2
    trails[0].length.should == 3
    trails[0][0].node_id.should == 3
    trails[0][1].node_id.should == 2
    trails[0][2].node_id.should == 1
    trails[1][0].node_id.should == 3
    trails[1][1].node_id.should == 4
    trails[1][2].node_id.should == 1
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
    initial_node = graph.nodes[3]
    terminal_node = graph.nodes[4]

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    GraphTesting.sorted_paths(trails).should == [[3,4]]

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 250, true)
    GraphTesting.sorted_paths(trails).should == []


    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[4]

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[4]
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[7]
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[8]
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[9]
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[9]
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, nil, true)
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
    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[4]
    graph.nodes[1].ends_of_kmers_of_node = 'A'*10
    graph.nodes[2].ends_of_kmers_of_node = 'A'*100
    graph.nodes[3].ends_of_kmers_of_node = 'A'*10
    graph.nodes[4].ends_of_kmers_of_node = 'A'*10
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new
    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 200, true)
    GraphTesting.sorted_paths(paths).should == [
      [1,2,3,4],
      [1,3,4],
    ].sort

    paths = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 100, true)
    (1..4).each do |node_id|
      graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
    end
    GraphTesting.sorted_paths(paths).should == [
      [1,3,4],
    ].sort

  end
end

# Dead code below I think
#describe 'Trail' do
#  it 'should give the right length' do
#    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
#
#    trail = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder::Trail.new
#    trail[0] = graph.nodes[1]
#    trail.nucleotide_length.should == 228+30
#    trail[1] = graph.nodes[2]
#    trail.nucleotide_length.should == 228+29+30
#    trail[2] = graph.nodes[4]
#    trail.nucleotide_length.should == 228+29+224+30
#  end
#end

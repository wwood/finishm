require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
class Util
  def self.revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

describe "AcyclicConnectionFinder" do
  it 'should calculate trails' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
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
    trails[0][1].node_id.should == 2
    trails[0][2].node_id.should == 4
    trails[1][0].node_id.should == 1
    trails[1][1].node_id.should == 3
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
    trails.length.should == 1
    trails[0].length.should == 2
    trails[0][0].node_id.should == 3
    trails[0][1].node_id.should == 4

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 250, true)
    trails.length.should == 0


    initial_node = graph.nodes[1]
    terminal_node = graph.nodes[4]

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 99999999, true)
    trails.length.should == 2
    trails[0].length.should == 3

    trails = cartographer.find_trails_between_nodes(graph, initial_node, terminal_node, 228+29+224+30+2, true)
    trails.length.should == 1
    trails[0].length.should == 3
    trails[0][0].node_id.should == 1
    trails[0][1].node_id.should == 2
    trails[0][2].node_id.should == 4
  end
end

describe 'Trail' do
  it 'should give the right length' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')

    trail = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder::Trail.new
    trail[0] = graph.nodes[1]
    trail.nucleotide_length.should == 228+30
    trail[1] = graph.nodes[2]
    trail.nucleotide_length.should == 228+29+30
    trail[2] = graph.nodes[4]
    trail.nucleotide_length.should == 228+29+224+30
  end
end

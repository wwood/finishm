require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "coverage based graph filter" do
  it 'should filter nodes by coverage' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')

    graph.nodes.length.should == 4
    cutoffer = Bio::AssemblyGraphAlgorithms::CoverageBasedGraphFilter.new

    deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, 0)
    deleted_nodes.length.should == 0
    deleted_arcs.length.should == 0
    graph.nodes.length.should == 4

    deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, 3.5)
    deleted_nodes.length.should == 2
    deleted_arcs.length.should == 4
    graph.nodes.length.should == 2

    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, 10)
    deleted_nodes.length.should == 4
    deleted_arcs.length.should == 4
    graph.nodes.length.should == 0
  end

  it 'should respect whitelisted sequences' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    graph.nodes.length.should == 4
    cutoffer = Bio::AssemblyGraphAlgorithms::CoverageBasedGraphFilter.new

    deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, 10, :whitelisted_sequences => [1])
    deleted_nodes.length.should == 1
    graph.nodes.length.should == 3

    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    deleted_nodes, deleted_arcs = cutoffer.remove_low_coverage_nodes(graph, 10, :whitelisted_sequences => [1,2])
    deleted_nodes.length.should == 0
    deleted_arcs.length.should == 0
    graph.nodes.length.should == 4
  end

  it 'should filter by connectivity' do
    graph = GraphTesting.emit([
      [1,2],
      [3,4]
    ])
    filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
    graph.nodes.length.should == 4
    filter.remove_unconnected_nodes(graph,[graph.nodes[1]])
    graph.nodes.length.should == 2
    graph.arcs.length.should == 1

    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4]
    ])
    filter.remove_unconnected_nodes(graph,[graph.nodes[1]])
    graph.nodes.length.should == 4
  end

  it 'should filter by leashed connectivity' do
    graph = GraphTesting.emit([
      [1,2],
      [3,4]
    ])
    filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
    graph.nodes.length.should == 4
    filter.remove_unconnected_nodes(graph,[graph.nodes[1]], :leash_length => 5)
    graph.nodes.length.should == 2 #or 1.meh.
    graph.arcs.length.should == 1
  end

  it 'should filter by leashed connectivity 2 nodes' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4]
    ])
    graph.nodes.length.should == 4
    filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
    filter.remove_unconnected_nodes(graph,[graph.nodes[1]], :leash_length => 15)
    graph.nodes.length.should == 3#or 2.meh.

  end

  it 'should filter by fake leashed connectivity' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,4]
    ])
    graph.nodes.length.should == 4
    filter = Bio::AssemblyGraphAlgorithms::ConnectivityBasedGraphFilter.new
    filter.remove_unconnected_nodes(graph,[graph.nodes[1]], :leash_length => 9999)
    graph.nodes.length.should == 4

  end
end

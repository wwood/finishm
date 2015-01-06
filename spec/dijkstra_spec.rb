require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "Dijkstra" do

  "SOME OF THIS IS LIKELY TESTED IN OTHER FILES SO GO LOOKING BEFORE NEW TESTS ARE WRITTEN"

  it 'should do hello world' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2]])
    GraphTesting.make_reads_paired(f)
    dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
    dijkstra.min_distances(f.graph, onode).should == {[1, :start_is_first]=>0, [2, :start_is_first]=>0}
  end

  it 'should do paired end jumps' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2],[4,2]])
    GraphTesting.make_reads_paired(f)
    dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)

    dijkstra.min_distances(f.graph, onode, :neighbour_finder => finder
      ).should == {[1, :start_is_first]=>0,
        [2, :start_is_first]=>0,
        [3, :end_is_first] => 80.0,
        [4, :end_is_first] => 80.0}
  end
end

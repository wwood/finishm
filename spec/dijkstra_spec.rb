require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "Dijkstra" do
  describe 'min distances' do

    "SOME OF THIS IS LIKELY TESTED IN OTHER FILES SO GO LOOKING BEFORE NEW TESTS ARE WRITTEN"

    it 'should do hello world' do
      f = GraphTesting.finishm_graph([[1,2],[3,4]])
      GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2]])
      GraphTesting.make_reads_paired(f)
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
      # won't find node 2 because the default finder is single-ended
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

    it 'should respect max-nodes' do
      f = GraphTesting.finishm_graph([[1,2],[2,3],[3,4]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
      dijkstra.min_distances(f.graph, onode).should == {
        [1, :start_is_first]=>0,
        [2, :start_is_first]=>0,
        [3, :start_is_first]=>10,
        [4, :start_is_first]=>20,
        }
      dijkstra.min_distances(f.graph, onode, :max_nodes => 2).should == {
        [1, :start_is_first]=>0,
        [2, :start_is_first]=>0,
        [3, :start_is_first]=>10,
        }
      dijkstra.min_distances(f.graph, onode, :max_nodes => 0).should == {
        [1, :start_is_first]=>0,
        }
      dijkstra.min_distances(f.graph, onode, :max_nodes => 3).should == {
        [1, :start_is_first]=>0,
        [2, :start_is_first]=>0,
        [3, :start_is_first]=>10,
        [4, :start_is_first]=>20,
        }
      dijkstra.min_distances(f.graph, onode, :max_nodes => 30).should == {
        [1, :start_is_first]=>0,
        [2, :start_is_first]=>0,
        [3, :start_is_first]=>10,
        [4, :start_is_first]=>20,
        }
    end

    it 'should max-nodes when ignoring directions' do
      f = GraphTesting.finishm_graph([[1,2],[2,3],[3,4]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)

      dijkstra.min_distances(f.graph, onode, {:max_nodes => 0, :ignore_directions => true}).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        }
      dijkstra.min_distances(f.graph, onode, {:max_nodes => 1, :ignore_directions => true}).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        [2, :start_is_first]=>0,
        [2, :end_is_first]=>0,
        }
      dijkstra.min_distances(f.graph, onode, {:max_nodes => 2, :ignore_directions => true}).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        [2, :start_is_first]=>0,
        [2, :end_is_first]=>0,
        [3, :start_is_first]=>10,
        [3, :end_is_first]=>10
        }
    end

    it 'should be ok with max-nodes when there is a tie' do
      f = GraphTesting.finishm_graph([[1,2],[1,3]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)

      dijkstra.min_distances(f.graph, onode, {:max_nodes => 0, :ignore_directions => true}).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        }
      dijkstra.min_distances(f.graph, onode, {:max_nodes => 1, :ignore_directions => true}).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        [2, :start_is_first]=>0,
        [2, :end_is_first]=>0,
        [3, :start_is_first]=>0,
        [3, :end_is_first]=>0,
        }

      f = GraphTesting.finishm_graph([[10,1],[1,2],[1,3]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[10], :start_is_first)
      dijkstra.min_distances(f.graph, onode, {:max_nodes => 2, :ignore_directions => true}).should == {
        [10, :start_is_first]=>0,
        [10, :end_is_first]=>0,
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        [2, :start_is_first]=>10,
        [2, :end_is_first]=>10,
        [3, :start_is_first]=>10,
        [3, :end_is_first]=>10,
        }
    end

    it 'should not get tripped up on negative distances' do
      f = GraphTesting.finishm_graph([[1,2],[3,4]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
      GraphTesting.add_reads_to_nodes(f, [[2,1],[3,2]])
      GraphTesting.make_reads_paired(f)
      GraphTesting.set_node_length(f.graph.nodes[2], 1000)
      finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
      dijkstra.min_distances(f.graph, onode, {
        :max_nodes => 1, :ignore_directions => true, :neighbour_finder => finder
        }).should == {
        [1, :start_is_first]=>0,
        [1, :end_is_first]=>0,
        [2, :start_is_first]=>0,
        [2, :end_is_first]=>0,
        }
    end
  end

  describe 'trail stats' do
    it 'should return nodes grouped by distance' do
      f = GraphTesting.finishm_graph([[1,2],[1,3],[2,3]])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
      dijkstra.trail_stats(f.graph, onode).should == {
        0=>{
          [1, :start_is_first]=>1,
          [2, :start_is_first]=>1,
          [3, :start_is_first]=>1,
          },
        10=>{
          [3, :start_is_first]=>1
          }
        }
    end

    it 'should track node positional multiplicities' do
      f = GraphTesting.finishm_graph([
        [1,3],
        [2,3],
        [3,4],
        [3,5],
        [4,6],
        [5,6],
        [6,7]
        ])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[3], :start_is_first)
      dijkstra.trail_stats(f.graph, onode).should == {
        0=>{
          [3, :start_is_first]=>1,
          [4, :start_is_first]=>1,
          [5, :start_is_first]=>1
          },
        10=>{
          [6, :start_is_first]=>2,
          },
        20=>{
          [7, :start_is_first]=>2,
          }
        }
    end

    it 'should deal with long divergent paths' do
      f = GraphTesting.finishm_graph([
        [1,2],
        [1,3],
        [2,4],
        [3,4],
        [4,5],
        [5,6]
        ])
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
      f.graph.nodes[3].ends_of_kmers_of_node = 'T'*100
      f.graph.nodes[5].ends_of_kmers_of_node = 'T'*100 #long paths
      dijkstra.trail_stats(f.graph, onode).should == {
        0=>{
          [1, :start_is_first]=>1,
          [2, :start_is_first]=>1,
          [3, :start_is_first]=>1,
          },
        10=>{
          [4, :start_is_first]=>1,
          },
        20=>{
          [5, :start_is_first]=>1,
          },
        100=>{
          [4, :start_is_first]=>1,
          },
        110=>{
          [5, :start_is_first]=>1,
          },
        120=>{
          [6, :start_is_first]=>1,
          },
        210=>{
          [6, :start_is_first]=>1,
          }
        }
    end

    it 'should ignore directions' do
      f = GraphTesting.finishm_graph([
        [1,3],
        [2,3],
        [3,4]
        ])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[3], :end_is_first)
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      dijkstra.trail_stats(f.graph, onode, {
        :ignore_directions => true, :max_nodes => 2
        }).should == {
        0=>{
          [3, :start_is_first]=>1,
          [3, :end_is_first]=>1,
          [1, :start_is_first]=>2,
          [1, :end_is_first]=>2,
          [2, :start_is_first]=>2,
          [2, :end_is_first]=>2
          }
        }
    end

    it 'should respect max nodes' do
      f = GraphTesting.finishm_graph([
        [1,3],
        [2,3],
        [3,4],
        [4,5],
        [4,6],
        [5,6]
        ])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[3], :start_is_first)
      dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      dijkstra.trail_stats(f.graph, onode, { :max_nodes => 1 }).should == {
        0=>{
          [3, :start_is_first]=>1,
          [4, :start_is_first]=>1
          }
        }

      dijkstra.trail_stats(f.graph, onode, { :max_nodes => 2 }).should == {
        0=>{
          [3, :start_is_first]=>1,
          [4, :start_is_first]=>1
          },
        10=>{
          [5, :start_is_first]=>1,
          [6, :start_is_first]=>1
          }
        }

      dijkstra.trail_stats(f.graph, onode, { :max_nodes => 3 }).should == {
        0=>{
          [3, :start_is_first]=>1,
          [4, :start_is_first]=>1
          },
        10=>{
          [5, :start_is_first]=>1,
          [6, :start_is_first]=>1
          }
        }

      dijkstra.trail_stats(f.graph, onode, { :max_nodes => 4 }).should == {
        0=>{
          [3, :start_is_first]=>1,
          [4, :start_is_first]=>1
          },
        10=>{
          [5, :start_is_first]=>1,
          [6, :start_is_first]=>1
          },
        20=>{
          [6, :start_is_first]=>1
          }
        }
    end
  end
end

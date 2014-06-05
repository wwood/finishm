require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "SingleCoherentWanderer" do
  describe 'wander' do
    it 'should calculate connections using a simple depth first search' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
      ])
      finishm_graph = Bio::FinishM::ProbedGraph.new
      finishm_graph.graph = graph
      finishm_graph.probe_nodes = [graph.nodes[1],graph.nodes[2]]
      finishm_graph.probe_node_directions = [true, false]

      cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new
      cartographer.wander(finishm_graph, 10000, 1, {}).should == {
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

      cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new
      cartographer.wander(finishm_graph, 10000, 1, {}).should == {
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

      cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new
      cartographer.wander(finishm_graph, 10000, 1, {}).should == {
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

      cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new
      cartographer.wander(finishm_graph, 20, 1, {}).should == {
        [0,1] => 20,
        #[0,3] => 40,
        [2,3] => 10,
      }
    end

    it 'should respect recoherence in the simplest case' do
      graph, initial_path, terminal_onode = GraphTesting.emit_ss([
        [1,2],

        [2,3],
        [3,4],
        [4,5],

        [2,6],
        [6,5],
      ], 1, 5)
      GraphTesting.add_noded_reads(graph,[
        [1,2,3,4,5],
        [2,6,5], #1,2,6 won't survive recoherence
        ])
      finishm_graph = Bio::FinishM::ProbedGraph.new
      finishm_graph.graph = graph
      finishm_graph.probe_nodes = [
        graph.nodes[1],graph.nodes[5],
      ]
      finishm_graph.probe_node_directions = [true, false]

      cartographer = Bio::AssemblyGraphAlgorithms::SingleCoherentWanderer.new
      cartographer.wander(finishm_graph, 100, 100, {}).should == {
        [0,1] => 40,
        } #with real recoherence have to take the longer way
      cartographer.wander(finishm_graph, 100, 1, {}).should == {
        [0,1] => 30,
        }#with no recoherence have to take the shorter way
    end
  end
end


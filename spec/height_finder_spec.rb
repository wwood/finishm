require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.array_of_sorted_nodes(nodes_array)
    nodes_array.collect do |nodes|
      nodes.collect{|n| n.node_id}.sort
    end
  end

  def self.array_of_sorted_nodes_in(nodes_array)
    to_return = {}
    nodes_array.flatten.each do |node|
      to_return[node.node_id] = node.nodes_in.collect{|n| n.node_id}.sort
    end
    return to_return
  end

  def self.array_of_sorted_nodes_out(nodes_array)
    to_return = {}
    nodes_array.flatten.each do |node|
      to_return[node.node_id] = node.nodes_out.collect{|n| n.node_id}.sort
    end
    return to_return
  end
end

describe "HeightFinder" do
  it 'should return maximum number of paths through simple graph' do
    f = GraphTesting.finishm_graph([[1,2],[1,3],[2,3]])
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], true) #forwards
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, [onode])
    GraphTesting.array_of_sorted_nodes(by_height).should == [
      [3],
      [2],
      [1]
      ]
    GraphTesting.array_of_sorted_nodes_out(by_height).should == {
      3=>[],
      2=>[3],
      1=>[2,3]
      }
    GraphTesting.array_of_sorted_nodes_in(by_height).should == {
      3=>[1,2],
      2=>[1],
      1=>[]
      }
    height_finder.max_paths_through(by_height).should == 2
    height_finder.min_paths_through(by_height).should == 1
    cycles.should == []
  end

  it 'should go in reverse' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [1,3],
      [2,3]
      ])
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[3], true # forwards
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, [onode], :reverse=>true)
    GraphTesting.array_of_sorted_nodes(by_height).should == [
      [1],
      [2],
      [3]
      ]
    GraphTesting.array_of_sorted_nodes_in(by_height).should == {
      1=>[2,3],
      2=>[3],
      3=>[]
      }
    GraphTesting.array_of_sorted_nodes_out(by_height).should == {
      1=>[],
      2=>[1],
      3=>[1,2]
      }
    height_finder.min_paths_through(by_height).should == 1
    height_finder.max_paths_through(by_height).should == 2
    cycles.should == []
  end

  it 'should respect range' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [1,3],
      [2,3],
      [3,4]
    ])
    range = [1, 3, 4].collect{|id| f.graph.nodes[id]}
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], true) # forwards
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, [onode], { :range=>range })
    GraphTesting.array_of_sorted_nodes(by_height).should == [
      [4],
      [3],
      [1]
      ]
    GraphTesting.array_of_sorted_nodes_in(by_height).should == {
      4=>[3],
      3=>[1],
      1=>[]
      }
    GraphTesting.array_of_sorted_nodes_out(by_height).should == {
      4=>[],
      3=>[4],
      1=>[3]
      }
  end

  it 'should cope with multiple initial nodes' do
    f = GraphTesting.finishm_graph([
      [1,3],
      [2,3],
      [3,4],
      [3,5],
      [4,6],
      [5,6],
      [6,7]
      ])
    onodes = [2, 7, 3, 1, 5].collect{|id| Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[id], true} #begin in random order
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, onodes)
    GraphTesting.array_of_sorted_nodes(by_height).should == [
      [7],
      [6],
      [4,5],
      [3],
      [1,2]
      ]
    GraphTesting.array_of_sorted_nodes_in(by_height).should == {
      7=>[6],
      6=>[4,5],
      5=>[3],
      4=>[3],
      3=>[1,2],
      2=>[],
      1=>[]
      }
    GraphTesting.array_of_sorted_nodes_out(by_height).should == {
      7=>[],
      6=>[7],
      5=>[6],
      4=>[6],
      3=>[4,5],
      2=>[3],
      1=>[3]
      }
    height_finder.max_paths_through(by_height).should == 4
    height_finder.min_paths_through(by_height).should == 2
    cycles.should == []
  end

  it 'should cope with multiple bubbling paths' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [1,3],
      [1,4],
      [2,8],
      [3,5],
      [4,5],
      [5,6],
      [5,7],
      [6,8],
      [7,8],
      ])
    onodes = [1, 2, 3].collect{|id| Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[id], true}
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, onodes)
    GraphTesting.array_of_sorted_nodes(by_height).should == [
      [8],
      [2,6,7],
      [5],
      [3,4],
      [1]
      ]
    GraphTesting.array_of_sorted_nodes_out(by_height).should == {
      8=>[],
      7=>[8],
      6=>[8],
      5=>[6,7],
      4=>[5],
      3=>[5],
      2=>[8],
      1=>[2,3,4]
      }
    GraphTesting.array_of_sorted_nodes_in(by_height).should == {
      8=>[2,6,7],
      7=>[5],
      6=>[5],
      5=>[3,4],
      4=>[1],
      3=>[1],
      2=>[1],
      1=>[]
      }
    height_finder.max_paths_through(by_height).should == 5
    height_finder.min_paths_through(by_height).should == 3
    cycles.should == []
  end

  it 'should report cycles and then ignore them' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [2,4],
      [3,2],
      [3,4]
      ])
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], true)
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse(f.graph, [onode])
    height_finder.max_paths_through(by_height).should == 2
    height_finder.min_paths_through(by_height).should == 1
    GraphTesting.array_of_sorted_nodes(cycles).should == [
      [2,3]
      ]
  end

  it 'should report nested cycles' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3],
      [3,2],
      [3,4],
      [4,2],
      [4,5]
      ])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(graph.nodes[1], true)
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      by_height, cycles = height_finder.traverse(graph, [onode])
      GraphTesting.array_of_sorted_nodes(cycles).should == [
        [2,3],
        [2,3,4]
        ]

  end

  describe 'find_oriented_edge_of_range' do
    it 'should find start and end points of a trail' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3]
        ])
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      start_and_end_onodes = height_finder.find_oriented_edge_of_range(graph)
      GraphTesting.array_of_sorted_nodes(start_and_end_onodes).sort.should == [
        [1],
        [3]
        ]
    end

    it 'should find multiple start and end points' do
      graph = GraphTesting.emit([
        [3,4],
        [3,5],
        [1,3],
        [2,3]
        ])
      start_and_end_onodes = Bio::AssemblyGraphAlgorithms::HeightFinder.new.find_oriented_edge_of_range(graph)
      GraphTesting.array_of_sorted_nodes(start_and_end_onodes).sort.should == [
        [1,2],
        [4,5]
        ]
    end

    it 'should respect range argument' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
        [2,4]
        ])
      nodes = [1,2,3].collect{|id| graph.nodes[id]}
      start_and_end_onodes = Bio::AssemblyGraphAlgorithms::HeightFinder.new.find_oriented_edge_of_range(graph, nodes)
      GraphTesting.array_of_sorted_nodes(start_and_end_onodes).sort.should == [
        [1],
        [3]
        ]
    end

    it 'should deal with cyclic paths' do
      graph = GraphTesting.emit([
        [1,3],
        [2,3],
        [3,4],
        [4,3],
        [4,5],
        [4,6]
        ])
      start_and_end_onodes = Bio::AssemblyGraphAlgorithms::HeightFinder.new.find_oriented_edge_of_range(graph)
      GraphTesting.array_of_sorted_nodes(start_and_end_onodes).sort.should == [
        [1,2],
        [5,6]
        ]
    end

    it 'should cope with a ring graph' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
        [3,1]
        ])
      start_and_end_nodes = Bio::AssemblyGraphAlgorithms::HeightFinder.new.find_oriented_edge_of_range(graph)
      start_and_end_nodes.should == [[],[]]
    end
  end
end

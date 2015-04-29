require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.ids_from_traversal_array(node_arr)
    node_arr.collect{|node| node.onode.node.node_id}.sort
  end

  def self.ids_from_list_of_arrays(list_of_node_arr)
    list_of_node_arr.collect{|node_arr| ids_from_traversal_array(node_arr)}
  end

  def self.ids_in_from_list_of_arrays(list_of_node_arr)
    to_return = {}
    list_of_node_arr.flatten.each do |node|
      to_return[node.onode.node.node_id] = ids_from_traversal_array(node.nodes_in)
    end
    return to_return
  end

  def self.ids_out_from_list_of_arrays(list_of_node_arr)
    to_return = {}
    list_of_node_arr.flatten.each do |node|
      to_return[node.onode.node.node_id] = ids_from_traversal_array(node.nodes_out)
    end
    return to_return
  end

  def self.ids_from_list_of_onode_arrays(list_of_onode_arr)
    list_of_onode_arr.collect{|onode_arr| onode_arr.collect{|onode| onode.node.node_id}.sort}
  end
end

describe "HeightFinder" do
  it 'should return maximum number of paths through simple graph' do
    f = GraphTesting.finishm_graph([[1,2],[1,3],[2,3]])
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], true) #forwards
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse f.graph, [onode]
    GraphTesting.ids_from_list_of_arrays(by_height).should == [
      [3],
      [2],
      [1]
      ]
    GraphTesting.ids_out_from_list_of_arrays(by_height).should == {
      3=>[],
      2=>[3],
      1=>[2,3]
      }
    GraphTesting.ids_in_from_list_of_arrays(by_height).should == {
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
    by_height, cycles = height_finder.traverse(f.graph, [onode], { :reverse=>true })
    GraphTesting.ids_from_list_of_arrays(by_height).should == [
      [1],
      [2],
      [3]
      ]
    GraphTesting.ids_in_from_list_of_arrays(by_height).should == {
      1=>[2,3],
      2=>[3],
      3=>[]
      }
    GraphTesting.ids_out_from_list_of_arrays(by_height).should == {
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
    range = [1, 3, 4].collect{|id| Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[id], true}
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse f.graph, range[0..0], { :range=>range }
    GraphTesting.ids_from_list_of_arrays(by_height).should == [
      [4],
      [3],
      [1]
      ]
    GraphTesting.ids_in_from_list_of_arrays(by_height).should == {
      4=>[3],
      3=>[1],
      1=>[]
      }
    GraphTesting.ids_out_from_list_of_arrays(by_height).should == {
      4=>[],
      3=>[4],
      1=>[3]
      }
  end

  it 'should cope with multiple roots' do
    f = GraphTesting.finishm_graph([
      [1,3],
      [2,3],
      [3,4],
      [3,5],
      [4,6],
      [5,6],
      [6,7]
      ]) # permute edge order to randomise start order
    onodes = [2, 7, 3, 1, 5].collect{|id| Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[id], true}
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height, cycles = height_finder.traverse f.graph, onodes
    GraphTesting.ids_from_list_of_arrays(by_height).should == [
      [7],
      [6],
      [4,5],
      [3],
      [1,2]
      ]
    GraphTesting.ids_in_from_list_of_arrays(by_height).should == {
      7=>[6],
      6=>[4,5],
      5=>[3],
      4=>[3],
      3=>[1,2],
      2=>[],
      1=>[]
      }
    GraphTesting.ids_out_from_list_of_arrays(by_height).should == {
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
    by_height, cycles = height_finder.traverse f.graph, onodes
    GraphTesting.ids_from_list_of_arrays(by_height).should == [
      [8],
      [2,6,7],
      [5],
      [3,4],
      [1]
      ]
    GraphTesting.ids_out_from_list_of_arrays(by_height).should == {
      8=>[],
      7=>[8],
      6=>[8],
      5=>[6,7],
      4=>[5],
      3=>[5],
      2=>[8],
      1=>[2,3,4]
      }
    GraphTesting.ids_in_from_list_of_arrays(by_height).should == {
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
    by_height, cycles = height_finder.traverse f.graph, [onode]
    height_finder.max_paths_through(by_height).should == 2
    height_finder.min_paths_through(by_height).should == 1
    GraphTesting.ids_from_list_of_onode_arrays(cycles).should == [
      [2,3]
      ]
  end

  it 'should report nested cycles' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,2],
      [3,4],
      [4,2],
      [4,5]
      ])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], true)
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      by_height, cycles = height_finder.traverse  f.graph, [onode]
      GraphTesting.ids_from_list_of_onode_arrays(cycles).should == [
        [2,3],
        [2,3,4]
        ]

  end

  describe 'find_bubbles' do
    it 'should report bubbles' do
      f = GraphTesting.finishm_graph([
        [1,2],
        [1,3],
        [2,4],
        [3,5],
        [3,6],
        [5,7],
        [6,7],
        [7,8],
        [7,9],
        [8,10],
        [9,10],
        [4,10],
        [10,11]
        ])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new f.graph.nodes[1], true
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      by_height, cycles = height_finder.traverse f.graph, [onode]
      bubbles = height_finder.find_bubbles(by_height)
      GraphTesting.ids_from_list_of_onode_arrays(bubbles).should == [
        [1,2,3,4,5,6,7,8,9,10],
        [3,5,6,7],
        [7,8,9,10] # this is pretty subtle bug
        ]
    end
  end

end

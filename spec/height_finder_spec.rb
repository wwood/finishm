require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.ids_from_node_array(node_arr)
    node_arr.collect{|node| node.onode.node.node_id}.sort
  end

  def self.ids_from_list_of_arrays(list_of_node_arr)
    list_of_node_arr.collect{|node_arr| ids_from_node_array(node_arr)}
  end

  def self.ids_in_from_list_of_arrays(list_of_node_arr)
    to_return = {}
    list_of_node_arr.flatten.each do |node|
      to_return[node.onode.node.node_id] = ids_from_node_array(node.nodes_in)
    end
    return to_return
  end

  def self.ids_out_from_list_of_arrays(list_of_node_arr)
    to_return = {}
    list_of_node_arr.flatten.each do |node|
      to_return[node.onode.node.node_id] = ids_from_node_array(node.nodes_out)
    end
    return to_return
  end
end

describe "HeightFinder" do
  it 'should return maximum number of paths through simple graph' do
    f = GraphTesting.finishm_graph([[1,2],[1,3],[2,3]])
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height = height_finder.traverse f.graph # forwards
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
  end

  it 'should cope with multiple roots' do
    f = GraphTesting.finishm_graph([
      [3,4],
      [4,6],
      [6,7],
      [3,5],
      [2,3],
      [1,3],
      [5,6]
      ]) # permute edge order to randomise start order
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height = height_finder.traverse f.graph
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
  end

  it 'should cope with multiple bubbling paths' do
    f = GraphTesting.finishm_graph([
      [3,5],
      [4,5],
      [5,6],
      [5,7],
      [6,8],
      [7,8],
      [1,2],
      [1,3],
      [1,4],
      [2,8]
      ])
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    by_height = height_finder.traverse f.graph
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
  end

  it 'should cope with cycles by arbitrarily ignoring linking paths' do
    # hmmm, this was definitely not specced to match the algorithm...
    f = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [2,4],
      [2,5],
      [3,2],
      [3,6],
      [3,7],
      [4,7],
      [5,7]
      ])
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      by_height = height_finder.traverse f.graph
      GraphTesting.ids_from_list_of_arrays(by_height).should == [
        [6,7],
        [4,5],
        [2],
        [1,3]
        ]
      height_finder.max_paths_through(by_height).should == 4 # fails
      height_finder.min_paths_through(by_height).should == 3 # fails

  end
end

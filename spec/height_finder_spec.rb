require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.ids_from_list_of_arrays(list_of_node_arr)
    list_of_node_arr.collect{|node_arr| node_arr.collect{|onode| onode.node.node_id}.sort}
  end

  def self.ids_from_hash_of_arrays(hash_of_node_arr)
    to_return = {}
    hash_of_node_arr.each_pair{|key, val| to_return[key[0]] = val.collect{|onode| onode.node.node_id}.sort}
    return to_return
  end

  def self.ids_from_hash_of_sets(hash_of_node_sets)
    to_return = {}
    hash_of_node_sets.each_pair{|key, val| to_return[key[0]] = val.collect{|settable| settable[0]}.sort}
    return to_return
  end
end

describe "HeightFinder" do
  it 'should return maximum number of paths through simple graph' do
    f = GraphTesting.finishm_graph([[1,2],[1,3],[2,3]])
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new f.graph
    height_finder.traverse
    GraphTesting.ids_from_list_of_arrays(height_finder.by_level).should == [
      [3],
      [2],
      [1]
      ]
    GraphTesting.ids_from_hash_of_arrays(height_finder.nodes_out).should == {
      2=>[3],
      1=>[2,3]
      }
    GraphTesting.ids_from_hash_of_sets(height_finder.nodes_in).should == {
      3=>[1,2],
      2=>[1]
      }
    height_finder.max_paths_through.should == 2
    height_finder.min_paths_through.should == 2
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
      ])
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    height_finder.traverse
    GraphTesting.ids_from_list_of_arrays(height_finder.by_level).should == [
      [7],
      [6],
      [4,5],
      [3],
      [1,2]
      ]
    GraphTesting.ids_from_hash_of_arrays(height_finder.nodes_out).should == {
      7=>[6],
      6=>[4,5],
      5=>[3],
      4=>[3],
      3=>[1,2]
      }
    GraphTesting.ids_from_hash_of_sets(height_finder.nodes_in).should == {
      6=>[7],
      5=>[6],
      4=>[6],
      3=>[4,5],
      2=>[3],
      1=>[3]
      }
    height_finder.max_paths_through.should == 4
    height_finder.min_paths_through.should == 2
  end

  it 'should cope with multiple bubbling paths' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [1,3],
      [1,4],
      [3,5],
      [4,5],
      [5,6],
      [5,7],
      [6,8],
      [7,8],
      [2,8]
      ])
    height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
    height_finder.traverse
    GraphTesting.ids_from_list_of_arrays(height_finder.by_level).should == [
      [8],
      [2,6,7],
      [5],
      [3,4],
      [1]
      ]
    GraphTesting.ids_from_hash_of_arrays(height_finder.nodes_out).should == {
      7=>[8],
      6=>[8],
      5=>[6,7],
      4=>[5],
      3=>[5],
      2=>[8],
      1=>[2,3,4]
      }
    GraphTesting.ids_from_hash_of_sets(height_finder.nodes_in).should == {
      8=>[2,6,7],
      7=>[5],
      6=>[5],
      5=>[3,4],
      4=>[1],
      3=>[1],
      2=>[1]
      }
    height_finder.max_paths_through.should == 5
    height_finder.min_paths_through.should == 3
  end

  it 'should cope with cycles by arbitrarily ignoring linking paths' do
    f = GraphTesting.finishm_graph([
      [1,2],
      [1,3],
      [1,4],
      [2,1],
      [2,5],
      [2,6],
      [3,5],
      [4,6]
      ])
      height_finder = Bio::AssemblyGraphAlgorithms::HeightFinder.new
      height_finder.traverse
      height_finder.max_paths_through.should == 4
      height_finder.min_paths_through.should == 3

  end
end

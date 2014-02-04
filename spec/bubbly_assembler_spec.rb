require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

class GraphTesting
  def self.metapath_to_array(metapath)
    to_return = []
    metapath.collect do |node_or_arr|
      if node_or_arr.kind_of?(Array)
        to_return.push node_or_arr.collect{|onode| onode.node_id}
      elsif node_or_arr.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        to_return.push node_or_arr.node_id
      else
        raise "Unknown metapath element: #{node_or_arr.inspect}"
      end
    end
    return to_return
  end
end


describe "Bubble Assembler" do
  it 'should assemble something easy' do
    graph, initial_path, terminal = GraphTesting.emit_ss([
      [1,2],
      [2,3],
    ], 1, 3)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new
    metapath = cartographer.assemble_from_node(graph, initial_path, nil)
    GraphTesting.metapath_to_array(metapath).should == [1,2,3]
  end

  it 'should handle a simple bubble' do
    graph, initial_path, terminal = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      [3,5],
    ], 1, 3)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new
    metapath = cartographer.assemble_from_node(graph, initial_path, nil)
    GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5]
  end
end




describe 'metapath' do
  it 'should act like an array' do
    metapath = Bio::AssemblyGraphAlgorithms::BubblyAssembler::MetaPath.new
    metapath.length.should == 0
    graph, initial_path, term = GraphTesting.emit_ss([[1,2]],1,2)
    onode = initial_path[0]
    metapath << onode
    metapath.length.should == 1
    metapath[0].should == onode
  end
end
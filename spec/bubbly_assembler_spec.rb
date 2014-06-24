require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

class GraphTesting
  def self.metapath_to_array(metapath)
    to_return = []
    metapath.collect do |node_or_arr|
      if node_or_arr.kind_of?(Array)
        to_return.push node_or_arr.collect{|otrail| otrail.length == 1 ? otrail[0].node_id : otrail.collect{|onode| onode.node_id}}
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
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
    metapath = cartographer.assemble_from_node(initial_path, nil)
    GraphTesting.metapath_to_array(metapath).should == [1,2,3]
  end

  it 'should handle a simple bubble' do
    graph, initial_path, terminal = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      [3,5],
    ], 1, 5)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
    metapath = cartographer.assemble_from_node(initial_path, nil)
    GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5]
  end

  it 'should handle 3 way bubbles that don\'t all converge on the same bubble' do
    graph, initial_path, terminal = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [3,4],
      [4,5],
      [2,6],
      [6,7],
      [7,5],
      [6,4],
      [7,8],
    ], 1, 8)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
    metapath = cartographer.assemble_from_node(initial_path, nil)
    GraphTesting.metapath_to_array(metapath).should == [1,2,[[3,4],[6,4],[6,7]],5,6]
  end

  it 'should deal with fluff in the middle of a bubble' do
    fail
  end

  it 'should be able to get out of the middle of a simple path' do
  fail
  end

  it 'should be able to get out of the middle of a path when it is bubbly' do
  fail
  end

  it 'should be able to assemble several contigs' do
  fail
  end

  it 'should be able to use recoherence to get around an inter-genome repeat' do
  fail
  end

  it 'should respect the leash length' do
  fail
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

  it 'should be able to output a reference contig' do
  fail
  end

  it 'should be able to iterate over variants from the reference contig' do
  fail
  end

  it 'should be able to write out a reference and VCF format file' do
  fail
  end

  it 'should reverse!' do
    fail
  end
end

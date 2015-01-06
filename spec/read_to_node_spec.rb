require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "ReadToNode" do
  it 'should work with postive and negative node IDs' do
    rtn = Bio::FinishM::ReadToNode.new(File.join(TEST_DATA_DIR, 'read_to_node/1_a_graph/ReadToNode.bin'))
    rtn[1].should == [1]
    rtn[2].should == [1]
    rtn[256].should == [7]
    rtn[6002].should == [18]
  end

  it 'should be able to deal with reads that are not contained in the graph at all' do
    rtn = Bio::FinishM::ReadToNode.new(File.join(TEST_DATA_DIR, 'read_to_node/2_no_read256_or_259/ReadToNode.bin'))

    # Originally: read, node
    # 255, -3
    # 256, -7 => removed
    # 257, -1
    # 258, 2
    # 259, 1 => removed
    # 260, -1
    rtn[255].should == [3]
    rtn[256].should == []
    rtn[257].should == [1]
    rtn[258].should == [2]
    rtn[259].should == []
    rtn[260].should == [1]
  end

  it 'should be able to deal with there being no last read' do
    rtn = Bio::FinishM::ReadToNode.new(File.join(TEST_DATA_DIR, 'read_to_node/3_no_last_read/ReadToNode.bin'))
    rtn[6002].should == []
    rtn[6001].should == [3]
  end
end
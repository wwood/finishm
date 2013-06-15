require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

class String
  def revcom
    Bio::Sequence::NA.new(self).reverse_complement.to_s.upcase
  end
end

describe "OrientedNodeTrail" do
  it "should be able to store a sequence of oriented nodes" do
    graph = Bio::Velvet::Graph.parse_from_file File.expand_path("#{TEST_DATA_DIR}/velvet_test_trails/Assem/LastGraph")
    trail = Bio::Velvet::Graph::OrientedNodeTrail.new

    trail.to_a.length.should == 0
    trail.add_node graph.nodes[1], :start_is_first
    trail.to_a.length.should == 1
    trail.to_a[0].node.should == graph.nodes[1]
    trail.to_a[0].first_side.should == :start_is_first

    trail.add_node graph.nodes[2], :start_is_first
    trail.add_node graph.nodes[4], :end_is_first
    trail.to_a.length.should == 3
    trail.to_a[2].node.should == graph.nodes[4]
    trail.to_a[2].first_side.should == :end_is_first

    expect {trail.add_node graph.nodes[4], :no_side}.to raise_error
  end

  it 'should get the sequence of no nodes' do
    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.sequence.should == ''
  end

  it 'should get the sequence of one node' do
    graph = Bio::Velvet::Graph.parse_from_file File.expand_path("#{TEST_DATA_DIR}/velvet_test_trails/Assem/LastGraph")

    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.add_node graph.nodes[1], :start_is_first
    trail.sequence.should == graph.nodes[1].sequence
    trail.sequence

    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.add_node graph.nodes[1], :end_is_first
    trail.sequence.should == graph.nodes[1].sequence.revcom
  end


  it 'should get the sequence of three nodes' do
    graph = Bio::Velvet::Graph.parse_from_file File.expand_path("#{TEST_DATA_DIR}/velvet_test_trails/Assem/LastGraph")

    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.add_node graph.nodes[1], :start_is_first
    trail.add_node graph.nodes[2], :start_is_first
    trail.add_node graph.nodes[4], :end_is_first
    exp = 'CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
ACTATGCTGGTATTTCACTTCCAGGTACAGG'.gsub(/\n/,'')
    trail.sequence.should == exp
  end

  it 'should not do sequence right when there is not enough info' do
    graph = Bio::Velvet::Graph.parse_from_file File.expand_path("#{TEST_DATA_DIR}/velvet_test_trails/Assem/LastGraph")

    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.add_node graph.nodes[2], :start_is_first
    expect {trail.sequence}.to raise_error
  end
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')

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

#  # Below are two things I was testing, which I think might be a problem with velvet, but not sure.
#  it 'should not have the bug I found any more full' do
#    # This graph is somewhat abbreviated and not a full and complete LastGraph file, but should be ok
#    lastgraph = <<EOF
#5     2805    43      1
#NODE    1      71      0       0       13651   13188
#ATAGATTATTTTTATTTTTCAGAGAATTTACAGAAAGATCAGTTAAAATCCAGAGCAAGAAAAGCATTGCA
#AAATTCTCTGAAAAATAAAAATAATCTATGCTGCTCGGTTTGACAAATTTTATCCCGTAAACTCCCTTTTT
#NODE    2      47      0       0       9848    9462
#GAAATTAAAAGAAGCTAAAGAGATTCAAAAATCGATCGATAACAAGA
#ATTTCTGCAATGCTTTTCTTGCTCTGGATTTTAACTGATCTTTCTGT
#NODE    3      11      0       0       2338    2246
#AACAGCTTCCA
#AGCTTCTTTTA
#NODE    4     2       0       0       396     378
#TA
#CA
#NODE    5     36      0       0       7023    6808
#TCATCAACAAGCTCAGCTTTGGATTTATCGAATTCT
#AGGAGTTTACGGGATAAAATTTGTCAAACCGAGCAG
#ARC     1      2      202
#EOF
#    lastgraph.gsub!(/ +/,"\t")
#    Tempfile.open('spec') do |f|
#      f.puts lastgraph
#      f.close
#
#      graph = Bio::Velvet::Graph.parse_from_file f.path
#      graph.nodes.length.should == 5
#      trail = Bio::Velvet::Graph::OrientedNodeTrail.new
#      trail.add_node graph.nodes[3], :end_is_first
#      trail.add_node graph.nodes[2], :end_is_first
#      trail.add_node graph.nodes[1], :end_is_first
#      trail.add_node graph.nodes[4], :start_is_first
#      trail.add_node graph.nodes[5], :start_is_first
#      trail.sequence.should == 'TGGAAGCTGTTTCTTGTTATCGATCGATTTTTGAATCTCTTTAGCTTCTTTTAATTTCTGCAATGCTTTTCTTGCTCTGGATTTTAACTGATCTTTCTGTAAATTCTCTGAAAAATAAAAATAATCTAT GCTGCTCGGTTTGACAAATTTTATCCCGTAAACTCCTTTTTTATCATCAACAAGCTCAGCTTTGGATTTATCGAATTCT'
#    end
#  end
#
#  it 'should not have the bug I found any more cut down version' do
#    # This graph is somewhat abbreviated and not a full and complete LastGraph file, but should be ok
#    lastgraph = <<EOF
#416     2705    43      1
#NODE    62      71      13651   13188   0       0
#ATAGATTATTTTTATTTTTCAGAGAATTTACAGAAAGATCAGTTAAAATCCAGAGCAAGAAAAGCATTGCA
#AAATTCTCTGAAAAATAAAAATAATCTATGCTGCTCGGTTTGACAAATTTTATCCCGTAAACTCCCTTTTT
#NODE    165     2       396     378     0       0
#TA
#CA
#ARC	-62	165	198
#EOF
#    lastgraph.gsub!(/ +/,"\t")
#    Tempfile.open('spec') do |f|
#      f.puts lastgraph
#      f.close
#
#      graph = Bio::Velvet::Graph.parse_from_file f.path
#      graph.nodes.length.should == 2
#      trail = Bio::Velvet::Graph::OrientedNodeTrail.new
#      trail.add_node graph.nodes[62], :end_is_first
#      trail.add_node graph.nodes[165], :start_is_first
#      trail.sequence.should == 'TGCAATGCTTTTCTTGCTCTGGATTTTAACTGATCTTTCTGTAAATTCTCTGAAAAATAAAAATAATCTAT GCTGCTCGGTTTGACAAATTTTATCCCGTAAACTCCTTTTTTA'
#    end
#  end

  it 'should give only 1 direction when entering a 2 node loop' do
    graph = GraphTesting.emit([
      [1,2],
      [2,4],
      [4,2],
    ])
    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    trail.add_node graph.nodes[1], :start_is_first
    trail.add_node graph.nodes[2], :start_is_first
    trail.neighbours_of_last_node(graph).collect{|n| [n.node.node_id, n.first_side]}.should == [[4,:start_is_first]]
  end
end

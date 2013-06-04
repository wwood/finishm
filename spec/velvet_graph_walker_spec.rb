require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "velvet graph walker" do
  it 'should calculate trail sequences' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_trails','Assem','LastGraph')
    hash_length = 31
    # Make sure we are on the same page
    graph.arcs.length.should eq(4)
    graph.nodes.length.should eq(4)
    nodes = graph.nodes

    walker = Bio::AssemblyGraphAlgorithms::LazyGraphWalker.new

    #    nodes.each do |node|
    #      puts '>'+node.node_id.to_s+'_fwd'
    #      puts node.ends_of_kmers_of_node
    #      puts '>'+node.node_id.to_s+'_rev'
    #      puts node.ends_of_kmers_of_twin_node
    #    end

    Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

    # sequenceChop.pl read1.fa 1 258
    exp = 'CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAGTAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAGATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATACGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATGGACGAGTTATATTTACTG'
    walker.trail_sequence(graph, [nodes[1]]).should eq(exp)

    # sequenceChop.pl read1.fa 1 287
    exp = 'CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAGTAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAGATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATACGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATGGACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAA'
    walker.trail_sequence(graph, [nodes[1],nodes[2]]).should eq(exp)

    lambda {walker.trail_sequence(graph, [nodes[1],nodes[2],nodes[3]]).should raise_error(Bio::AssemblyGraphAlgorithms::GraphWalkingException)} #Nodes 2 and 3 don't share an arc

    # read1.fa all
    exp = 'CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAGTAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAGATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATACGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATGGACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTCCTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATGATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAAGTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAAACTATGCTGGTATTTCACTTCCAGGTACAGG'
    walker.trail_sequence(graph, [nodes[1],nodes[2],nodes[4]]).should eq(exp)

    # sequenceChop.pl read1.fa 1 287 |reverseFasta.pl
    exp = 'TTTTATAAAGTAATCTCCTTCTTTTAAACCAGTAAATATAACTCGTCCATTTTTATCAGTTACACCCTTTCCTTTTAATAAAACCACATTTCCAGTAGAATCATACGTATATTTACCAATTACATTACCATTTTTATCCCTAACAGAAAAAGCTGCGCCTGCAAGATCTATTGAAATATTTTCTGAATCTACTTTTTTAACTCCGAATCCCCATGTATAAGTTGTTACTTTATCTTCTAAAACTTTATAGTTTGATTCTAAATCGTGATCTTTGGTAGAGATAAGTG'
    #WRONG sequenceChop.pl spec/data/velvet_test_trails/read1.fa 1 288 |reverseFasta.pl
    #exp = 'CTTTTATAAAGTAATCTCCTTCTTTTAAACCAGTAAATATAACTCGTCCATTTTTATCAGTTACACCCTTTCCTTTTAATAAAACCACATTTCCAGTAGAATCATACGTATATTTACCAATTACATTACCATTTTTATCCCTAACAGAAAAAGCTGCGCCTGCAAGATCTATTGAAATATTTTCTGAATCTACTTTTTTAACTCCGAATCCCCATGTATAAGTTGTTACTTTATCTTCTAAAACTTTATAGTTTGATTCTAAATCGTGATCTTTGGTAGAGATAAGTG'
    walker.trail_sequence(graph, [nodes[2],nodes[1]]).should eq(exp)

    lambda {walker.trail_sequence(graph, [nodes[2],nodes[4],nodes[3]]).should raise_error(Bio::AssemblyGraphAlgorithms::GraphWalkingException)} #share arcs, but in the wrong direction (2 and 3 both go to the start of 4)

    walker.trail_sequence([]).should eq('')
  end
end

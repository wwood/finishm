require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')


describe "AllOrfs" do

  describe 'should find a hello world trail' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3]
      ])
    initial_path = GraphTesting.make_onodes(graph, %w(1s))
    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, initial_path)

    p problems
    orfs = orfer.find_orfs_from_problems(problems)

  end

  describe 'search_for_codons' do
    it 'should look for forward start codons and twin stop codons ending in first node of trail' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'AAATGAAAAA' # start codon 'ATG'
      graph.nodes[1].ends_of_kmers_of_twin_node = 'TTTTTAACTT' # stop codon 'TAA'
      orfer.search_for_codons(otrails[0], false).should == [
        [5],
        [7]
        ]
    end

    it 'should look for forward stop codons and twin start codons ending in last node of trail' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[2].ends_of_kmers_of_node = 'AAATAATAAA' # stop codon 'TAA'
      graph.nodes[2].ends_of_kmers_of_twin_node = 'TTGATGTTTT' # start codon 'ATG', stop codon 'TGA'
      orfer.search_for_codons(otrails[0], true).should == [
        [6,9],
        [6]
        ]
    end
  end

  describe 'get_overlap_sequences' do
    it 'should traverse trail to get enough sequence' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[2].ends_of_kmers_of_node = 'G'*10 # Second node is sequence of G's
      graph.nodes[2].ends_of_kmers_of_twin_node = 'C'*10
      orfer.get_overlap_sequences(otrails[0], 3).should == [
        'AAGG',
        'CCTT'
        ]
    end

    it 'should look across multiple nodes if necessary' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2,3]])
      graph.nodes[2].ends_of_kmers_of_node = 'G' # Second node is single G
      graph.nodes[2].ends_of_kmers_of_twin_node = 'C'
      orfer.get_overlap_sequences(otrails[0], 4).should == [
        'AAAGAA',
        'TTCTTT'
        ]
    end

    it 'should handle a short initial node' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'C' # First node is single C
      graph.nodes[1].ends_of_kmers_of_twin_node = 'G'
      orfer.get_overlap_sequences(otrails[0], 4).should == [
        'CAAA',
        'TTTG'
        ]
    end

    it 'should be able to work back from end of trail' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph, otrails = GraphTesting.emit_otrails([[1,2,3]])
      graph.nodes[2].ends_of_kmers_of_node = 'C'
      graph.nodes[2].ends_of_kmers_of_twin_node = 'G'
      orfer.get_overlap_sequences(otrails[0], 5, true).should == [
        'AAACAAAA',
        'TTTTGTTT'
        ]
    end
  end

  describe 'get_sequences' do
    it 'should get forward and twin sequences of a node alone' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      graph = GraphTesting.emit([[1,2]])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new graph.nodes[1], true

      orfer.get_sequences(onode).should == [
        'A'*10,
        'T'*10
        ]
    end
  end

  describe 'word_search' do
    it 'should report end position of words within strings' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.word_search('Hello', ['lo'], 2).should == {
        'lo' => [5]
        }
      orfer.word_search('Agitate, infiltrate', ['ate', 'lat'], 3).should == {
        'ate' => [7,19]
        }
    end
  end

  describe 'orfs_from_start_stop_indices' do
    it 'should work when there are no orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([],[],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.final_stop_positions.should == []
    end

    it 'should work for one orf' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0],[6],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[0,6]]
      res.final_stop_positions.should == []
    end

    it 'should work for one orf in 2 frames' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0,2],[6,11],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[0,6],[2,11]]
      res.final_stop_positions.should == []
    end

    it 'should work for one orf and a leftover' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([6],[0,9],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[6,9]]
      res.final_stop_positions.should == [0]
    end

    it 'should work with a start orf' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([],[7],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.final_stop_positions.should == [7]
    end

    it 'should work with 3 orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[11,17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[8,11],[14,17],[20,44]]
      res.final_stop_positions.should == []
    end

    it 'should work with an internal start codon' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[8,17],[20,44]]
      res.final_stop_positions.should == []
    end
  end
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')

class GraphTesting
  def self.sorted_path_results(paths, forwards=true)
    paths.sort_by{|path| path.collect{|n| n.node_id}}.collect do |path|
      forwards ? path.fwd_orfs_result : path.twin_orfs_result
    end
  end
end

describe "AllOrfs" do

  it 'should find a hello world ORF' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'AAATGGAAAA' #start codon 'ATG'
    graph.nodes[3].ends_of_kmers_of_node = 'AAAAAATAAA' #stop codon 'TAA'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))


    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, initial_path)
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    res.start_stop_pairs.should == [
      [5,29]
      ]
    res.initial_stop_positions.should == []
    res.final_start_positions.should == []
  end

  it 'should find ORFs over a bubble' do
    graph, initial_path, = GraphTesting.emit_ss([
      [1,2],
      [1,3],
      [2,4],
      [3,4]
      ], 1, 4)
    graph.nodes[1].ends_of_kmers_of_node = 'AAATGGAAAA' # start 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'C'
    graph.nodes[3].ends_of_kmers_of_node = 'A'
    graph.nodes[4].ends_of_kmers_of_node = 'AATTAAAAAA' # stop 'TAA'

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, initial_path)
    #pp problems
    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,4],
      [1,3,4]
      ]
    res = GraphTesting.sorted_path_results(paths.trails, true) # forward direction
    res.collect{|result| result.start_stop_pairs.sort}.should == [
      [[5,17]],
      [[5,17]]
      ]
    res.collect{|result| result.final_start_positions}.should == [[],[]]
    res.collect{|result| result.initial_stop_positions}.should == [[],[]]
  end

  it 'should respect phase along each trail' do
    graph, initial_path, = GraphTesting.emit_ss([
      [1,2],
      [1,3],
      [2,4],
      [3,4]
      ], 1, 4)
    graph.nodes[1].ends_of_kmers_of_node = 'AAATGGAAAA' # start 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'C'
    graph.nodes[3].ends_of_kmers_of_node = 'AAA'
    graph.nodes[4].ends_of_kmers_of_node = 'AATTAAAAAA' # stop 'TAA'

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, initial_path)
    #pp problems
    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,4],
      [1,3,4]
      ]
    res = GraphTesting.sorted_path_results(paths.trails, true) # forwards direction
    res.collect{|result| result.start_stop_pairs.sort}.should == [
      [[5,17]],
      []
      ]
    res.collect{|result| result.final_start_positions}.should == [
      [],
      [5]
      ]
    res.collect{|result| result.initial_stop_positions}.should == [
      [],
      [19]
      ]
  end

  it 'should respect leash length' do
    graph, initial_path, = GraphTesting.emit_ss([
      [1,2],
      [2,3]
      ], 1, 3)
    graph.nodes[1].ends_of_kmers_of_node = 'AAATGGAAAA' # start 'ATG'
    graph.nodes[3].ends_of_kmers_of_node = 'AAAAAATAGA' # stop 'TAG'

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, initial_path, :leash_length => 20)
    paths = orfer.find_orfs_from_problems(problems)
    GraphTesting.sorted_paths(paths.trails).should == []

    problems = orfer.find_all_problems(graph, initial_path, :leash_length => 40)
    paths = orfer.find_orfs_from_problems(problems)
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    res.start_stop_pairs.should == [[5,29]]
    res.initial_stop_positions.should == []
    res.final_start_positions.should == []
  end

  it 'should respect terminal nodes' do
    fail '#todo'
  end

  it 'should respect minimum orf length' do
    fail '#todo'
  end

  it 'should respect max gapfill paths' do
    fail '#todo'
  end

  it 'should respect max cycles' do
    fail '#todo'
  end

  it 'should respect max explore nodes' do
    fail '#todo'
  end

  describe 'search_for_codons' do
    it 'should report end positions for codons starting in first node of trail' do
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'AAATGAAAAA' # start codon 'ATG', stop codon 'TGA'
      graph.nodes[1].ends_of_kmers_of_twin_node = 'TTTTTAACTT' # stop codon 'TAA'

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      fwd_res, twin_res = orfer.search_for_codons(otrails[0])
      fwd_res.start_positions.should == [5]
      fwd_res.stop_positions.should == [6]
      twin_res.start_positions.should == []
      twin_res.stop_positions.should == [7]
    end

    it 'should work on single-node trail' do
      graph, otrails = GraphTesting.emit([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'AAATAATAAA' # stop codon 'TAA'
      graph.nodes[1].ends_of_kmers_of_twin_node = 'TTGATGTTTT' # start codon 'ATG', stop codon 'TGA'
      otrail = Bio::Velvet::Graph::OrientedNodeTrail.new
      otrail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      fwd_res, twin_res = orfer.search_for_codons(otrail)
      fwd_res.start_positions.should == []
      fwd_res.stop_positions.should == [6,9]
      twin_res.start_positions.should == [6]
      twin_res.stop_positions.should == [4]
    end
  end

  describe 'get_overlap_sequences' do
    it 'should traverse trail to get enough sequence' do
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[2].ends_of_kmers_of_node = 'G'*10 # Second node is sequence of G's
      graph.nodes[2].ends_of_kmers_of_twin_node = 'C'*10

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.get_overlap_sequences(otrails[0], 3).should == [
        'AAGG',
        'CCTT'
        ]
    end

    it 'should look across multiple nodes if necessary' do
      graph, otrails = GraphTesting.emit_otrails([[1,2,3]])
      graph.nodes[2].ends_of_kmers_of_node = 'G' # Second node is single G
      graph.nodes[2].ends_of_kmers_of_twin_node = 'C'

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.get_overlap_sequences(otrails[0], 4).should == [
        'AAAGAA',
        'TTCTTT'
        ]
    end

    it 'should handle a short initial node' do
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'C' # First node is single C
      graph.nodes[1].ends_of_kmers_of_twin_node = 'G'

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.get_overlap_sequences(otrails[0], 4).should == [
        'CAAA',
        'TTTG'
        ]
    end

    it 'should be able to work back from end of trail' do
      graph, otrails = GraphTesting.emit_otrails([[1,2,3]])
      graph.nodes[2].ends_of_kmers_of_node = 'C'
      graph.nodes[2].ends_of_kmers_of_twin_node = 'G'

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.get_overlap_sequences(otrails[0], 5, true).should == [
        'AAACAAAA',
        'TTTTGTTT'
        ]
    end
  end

  describe 'get_sequences' do
    it 'should get forward and twin sequences of a node alone' do
      graph = GraphTesting.emit([[1,2]])
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new graph.nodes[1], true

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
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
      res.initial_stop_positions.should == []
      res.final_start_positions.should == []
    end

    it 'should work for one orf' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0],[6],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[0,6]]
      res.initial_stop_positions.should == []
      res.final_start_positions.should == []
    end

    it 'should work for one orf in 2 frames' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0,2],[6,11],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[0,6],[2,11]]
      res.initial_stop_positions.should == []
      res.final_start_positions.should == []
    end

    it 'should work for one orf and a leftover' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([6],[0,9],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[6,9]]
      res.initial_stop_positions.should == [0]
      res.final_start_positions.should == []
    end

    it 'should work with unclosed orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([],[7],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.initial_stop_positions.should == [7]
      res.final_start_positions.should == []

      res = orfer.orfs_from_start_stop_indices([7],[],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.initial_stop_positions.should == []
      res.final_start_positions.should == [7]
    end

    it 'should work with 3 orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[11,17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[8,11],[14,17],[20,44]]
      res.initial_stop_positions.should == []
      res.final_start_positions.should == []
    end

    it 'should work with an internal start codon' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[8,17],[20,44]]
      res.initial_stop_positions.should == []
      res.final_start_positions.should == []
    end

    it 'should find first stop codon in a frame before an orf and first start after' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([7,13,19],[1,4,10],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[7,10]]
      res.initial_stop_positions.should == [1]
      res.final_start_positions.should == [13]
    end
  end
end

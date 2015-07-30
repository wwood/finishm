require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')

class GraphTesting
  def self.sorted_path_results(paths, forwards=true)
    paths.sort_by{|path| path.collect{|n| n.node_id}}.collect do |path|
      forwards ? path.fwd_orfs_result : path.twin_orfs_result
    end
  end

  def self.markers(start_positions, stop_positions)
    [start_positions, stop_positions].collect do |positions|
      positions.collect do |pos|
        marker = Bio::AssemblyGraphAlgorithms::AllOrfsFinder::Marker.new
        marker.position_in_trail = pos
        marker
      end
    end
  end

  def self.sorted_marker_pair_positions(pair_array)
    pair_array.collect do |pair|
      pair.collect{|m| m.position_in_trail}
    end.sort
  end

  def self.marker_positions(markers)
    markers.collect{|m| m.position_in_trail}
  end

  def self.sorted_marker_pair_node_positions(pair_array)
    pair_array.sort_by do |pair|
      pair.collect{|m| m.position_in_trail}
    end.collect do |pair|
      pair.collect{|m| [m.node.node_id, m.position_in_node]}
    end
  end

  def self.marker_node_positions(markers)
    markers.collect do |m|
      [m.node.node_id, m.position_in_node]
    end
  end
end

describe "AllOrfs" do

  it 'should find a hello world ORF' do
    graph, = GraphTesting.emit_otrails([
      [1,2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' #stop codon 'TAA', start codon 'ATG'
    graph.nodes[3].ends_of_kmers_of_node = 'AAAAAAATAA' #stop codon 'TAA'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [
      [6,30]
      ]
    GraphTesting.sorted_marker_pair_node_positions(res.start_stop_pairs).should == [
      [[1,6],[3,10]]
      ]
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [3]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[1,3]]
    res.final_start_markers.should == []
  end

  it 'should find a hello world ORF in twin direction' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_twin_node = 'TTTAGTTTTT' # stop codon 'TAG'
    graph.nodes[2].ends_of_kmers_of_twin_node = 'TAAATGTTTT' # stop codon 'TAA', start codon 'ATG'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].twin_orfs_result
    GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [
      [16,25]
      ]
    GraphTesting.sorted_marker_pair_node_positions(res.start_stop_pairs).should == [
      [[2,6],[1,5]]
      ]
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [13]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[2,3]]
    res.final_start_markers.should == []
  end

  it 'should find ORFs over a bubble' do
    graph = GraphTesting.emit([
      [1,2],
      [1,3],
      [2,4],
      [3,4]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' # stop codon 'TAA', start 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'C'
    graph.nodes[3].ends_of_kmers_of_node = 'A'
    graph.nodes[4].ends_of_kmers_of_node = 'AAATTAAAAA' # stop 'TAA'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems
    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,4],
      [1,3,4]
      ]
    res = GraphTesting.sorted_path_results(paths.trails, true) # forward direction
    res.collect{|result| GraphTesting.sorted_marker_pair_positions(result.start_stop_pairs)}.should == [
      [[6,18]],
      [[6,18]]
      ]
    res.collect{|result| result.final_start_markers}.should == [[],[]]
    res.collect{|result| result.initial_start_markers}.should == [[],[]]
    res.collect{|result| GraphTesting.marker_positions(result.initial_stop_markers)}.should == [
      [3],
      [3]
      ]
    res.collect{|result| GraphTesting.marker_node_positions(result.initial_stop_markers)}.should == [
      [[1,3]],
      [[1,3]]
      ]
  end

  it 'should respect phase along each trail' do
    graph = GraphTesting.emit([
      [1,2],
      [1,3],
      [2,4],
      [3,4]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' # stop 'TAA', start 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'C'
    graph.nodes[3].ends_of_kmers_of_node = 'AAA'
    graph.nodes[4].ends_of_kmers_of_node = 'AAATTAGAAA' # stop 'TAG'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems
    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,4],
      [1,3,4]
      ]
    res = GraphTesting.sorted_path_results(paths.trails, true) # forwards direction
    res.collect{|result| GraphTesting.sorted_marker_pair_positions(result.start_stop_pairs)}.should == [
      [[6,18]],
      []
      ]
    GraphTesting.sorted_marker_pair_node_positions(res[0].start_stop_pairs).should == [
      [[1,6],[4,7]]
      ]
    res.collect{|result| GraphTesting.marker_positions(result.final_start_markers)}.should == [
      [],
      [6]
      ]
    GraphTesting.marker_node_positions(res[1].final_start_markers).should == [
      [1,6]
      ]
    res.collect{|result| GraphTesting.marker_positions(result.initial_stop_markers)}.should ==[
      [3],
      [3,20]
      ]
    res.collect{|result| GraphTesting.marker_node_positions(result.initial_stop_markers)}.should == [
      [[1,3]],
      [[1,3],[4,7]]
      ]
    res.collect{|result| result.initial_start_markers}.should == [[],[]]
  end

  it 'should find two same-phase orfs along a trail' do
    graph, = GraphTesting.emit_otrails([
      [1,2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' #stop codon 'TAA', start codon 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'AATAAATGGA' #stop codon 'TAA', start codon 'ATG'
    graph.nodes[3].ends_of_kmers_of_node = 'AAAAAAATAA' #stop codon 'TAA'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [
      [6,15],
      [18,30]
      ]
    GraphTesting.sorted_marker_pair_node_positions(res.start_stop_pairs).should == [
      [[1,6],[2,5]],
      [[2,8], [3,10]]
      ]
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [3]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[1,3]]
    res.final_start_markers.should == []
  end

  it 'should end orfs at first stop codon in forward direction' do
    graph, = GraphTesting.emit_otrails([
      [1,2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' #stop codon 'TAA', start codon 'ATG'
    graph.nodes[2].ends_of_kmers_of_node = 'AATAAAAAGA' #stop codon 'TAA'
    graph.nodes[3].ends_of_kmers_of_node = 'AAAAAAATAA' #stop codon 'TAA'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [
      [6,15]
      ]
    GraphTesting.sorted_marker_pair_node_positions(res.start_stop_pairs).should == [
      [[1,6],[2,5]]
      ]
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [3]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[1,3]]
    res.final_start_markers.should == []
  end

  it 'should end orfs at first stop codon in twin direction' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_twin_node = 'TTAGTTTTTT' # stop codon 'TAG'
    graph.nodes[2].ends_of_kmers_of_twin_node = 'TTTAGTTTTT' # stop codon 'TAG'
    graph.nodes[3].ends_of_kmers_of_twin_node = 'TAAATGTTTT' # stop codon 'TAA', start codon 'ATG'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].twin_orfs_result
    GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [
      [6,15]
      ]
    GraphTesting.sorted_marker_pair_node_positions(res.start_stop_pairs).should == [
      [[3,6],[2,5]]
      ]
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [3]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[3,3]]
    res.final_start_markers.should == []
  end

  it 'should return the first initial stop codon in forward direction' do
    graph = GraphTesting.emit([
      [1,2],
      [2,3]
      ])
    graph.nodes[1].ends_of_kmers_of_node = 'AAATAGAAAA' # stop codon 'TAG'
    graph.nodes[2].ends_of_kmers_of_node = 'AATAGAAAAA' # stop codon 'TAG'
    initial_path = GraphTesting.make_onodes(graph, %w(1s))

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    problems = orfer.find_all_problems(graph, [initial_path])
    #pp problems

    paths = orfer.find_orfs_from_problems(problems)
    #pp paths
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3]
      ]
    res = paths.trails[0].fwd_orfs_result
    res.start_stop_pairs.should == []
    res.initial_start_markers.should == []
    GraphTesting.marker_positions(res.initial_stop_markers).should == [6]
    GraphTesting.marker_node_positions(res.initial_stop_markers).should == [[1,6]]
    res.final_start_markers.should == []
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

  describe 'search_for_codons' do
    it 'should report end positions for codons starting in first node of trail' do
      graph, otrails = GraphTesting.emit_otrails([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'AAATGAAAAA' # start codon 'ATG', stop codon 'TGA'
      graph.nodes[1].ends_of_kmers_of_twin_node = 'TTTTTAACTT' # stop codon 'TAA'

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      fwd_res, twin_res = orfer.search_for_codons(otrails[0])
      fwd_res.start_markers.collect{|m| m.position_in_node}.should == [5]
      fwd_res.stop_markers.collect{|m| m.position_in_node}.should == [6]
      fwd_res.start_markers.all?{|n| n.node.node_id == 1}.should == true
      fwd_res.stop_markers.all?{|n| n.node.node_id == 1}.should == true
      twin_res.start_markers.should == []
      twin_res.stop_markers.collect{|m| m.position_in_node}.should == [7]
      twin_res.start_markers.all?{|n| n.node.node_id == 1}.should == true
      twin_res.stop_markers.all?{|n| n.node.node_id == 1}.should == true
    end

    it 'should work on single-node trail' do
      graph, otrails = GraphTesting.emit([[1,2]])
      graph.nodes[1].ends_of_kmers_of_node = 'AAATAATAAA' # stop codon 'TAA'
      graph.nodes[1].ends_of_kmers_of_twin_node = 'TTGATGTTTT' # start codon 'ATG', stop codon 'TGA'
      otrail = Bio::Velvet::Graph::OrientedNodeTrail.new
      otrail.add_node graph.nodes[1], Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      fwd_res, twin_res = orfer.search_for_codons(otrail)
      fwd_res.start_markers.should == []
      fwd_res.stop_markers.collect{|m| m.position_in_node}.should == [6,9]
      twin_res.start_markers.collect{|m| m.position_in_node}.should == [6]
      twin_res.stop_markers.collect{|m| m.position_in_node}.should == [4]
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

  describe 'orfs_from_start_stop_markers' do
    it 'should work when there are no orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_markers([],[],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.initial_start_markers.should == []
      res.initial_stop_markers.should == []
      res.final_start_markers.should == []
    end

    it 'should skip an orf before first stop' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers [0],[6]
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs == []
      GraphTesting.marker_positions(res.initial_start_markers).should == [0]
      GraphTesting.marker_positions(res.initial_stop_markers).should == [6]
      res.final_start_markers.should == []
    end

    it 'should find an orf after a stop' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers [6], [0,9]
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [[6,9]]
      res.initial_start_markers.should == []
      GraphTesting.marker_positions(res.initial_stop_markers).should == [0]
      res.final_start_markers.should == []
    end

    it 'should work for one orf in 2 frames' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers([3,5],[0,2,6,11])
      res = orfer.orfs_from_start_stop_markers(start, stop, 0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [[3,6],[5,11]]
      res.initial_start_markers.should == []
      GraphTesting.marker_positions(res.initial_stop_markers).should == [0,2]
      res.final_start_markers.should == []
    end

    it 'should work with unclosed orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers([],[7])
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      res.initial_start_markers.should == []
      GraphTesting.marker_positions(res.initial_stop_markers).should == [7]
      res.final_start_markers.should == []

      start, stop = GraphTesting.markers([7],[])
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == []
      GraphTesting.marker_positions(res.initial_start_markers).should == [7]
      res.initial_stop_markers.should == []
      res.final_start_markers.should == []
    end

    it 'should work with 3 orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers([8,14,20],[2,11,17,44])
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [[8,11],[14,17],[20,44]]
      res.initial_start_markers.should == []
      GraphTesting.marker_positions(res.initial_stop_markers).should == [2]
      res.final_start_markers.should == []
    end

    it 'should work with an internal start codon' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers([8,14,20],[5,17,44])
      res = orfer.orfs_from_start_stop_markers(start,stop,0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [[8,17],[20,44]]
      res.initial_start_markers.should == []
      GraphTesting.marker_positions(res.initial_stop_markers).should == [5]
      res.final_start_markers.should == []
    end

    it 'should find first stop codon in a frame before an orf and first start after' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      start, stop = GraphTesting.markers([7,13,19],[1,4,10])
      res = orfer.orfs_from_start_stop_markers(start, stop, 0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      GraphTesting.sorted_marker_pair_positions(res.start_stop_pairs).should == [[7,10]]
      GraphTesting.marker_positions(res.initial_stop_markers).should == [1]
      GraphTesting.marker_positions(res.final_start_markers).should == [13]
    end
  end

  describe 'orf_sequences_from_trails' do
    it 'should return orf sequences for a hello world orf' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3]
        ])
      graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' #stop codon 'TAA', start codon 'ATG'
      graph.nodes[3].ends_of_kmers_of_node = 'AAAAAAATAA' #stop codon 'TAA'
      initial_path = GraphTesting.make_onodes(graph, %w(1s))

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      problems = orfer.find_all_problems(graph, [initial_path])
      #pp problems

      paths = orfer.find_orfs_from_problems(problems)
      #pp paths
      orfer.orf_sequences_from_trails(paths.trails).should == [
        ['(1s:6),2s,(3s:10)', 'ATGGAAAAAAAAAAAAAAAAAAAATAA'],
        [',(1s:3)', 'TAA'],
        ['1s,2s,3s', 'T'*30],
        ['1s,2s,3s', 'T'*27],
        ['1s,2s,3s', 'T'*27]
        ]
    end

    it 'should respect minimum orf length' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3]
        ])
      graph.nodes[1].ends_of_kmers_of_node = 'TAAATGGAAA' #stop codon 'TAA', start codon 'ATG'
      graph.nodes[3].ends_of_kmers_of_node = 'AAAAAAATAA' #stop codon 'TAA'
      initial_path = GraphTesting.make_onodes(graph, %w(1s))

      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      problems = orfer.find_all_problems(graph, [initial_path])

      paths  =  orfer.find_orfs_from_problems(problems, :min_orf_length => 30)
      orfer.orf_sequences_from_trails(paths.trails, 30).should == [
        ['1s,2s,3s', 'T'*30]
        ]

      paths  =  orfer.find_orfs_from_problems(problems, :min_orf_length => 20)
      orfer.orf_sequences_from_trails(paths.trails, 20).should == [
        ['(1s:6),2s,(3s:10)', 'ATGGAAAAAAAAAAAAAAAAAAAATAA'],
        ['1s,2s,3s', 'T'*30],
        ['1s,2s,3s', 'T'*27],
        ['1s,2s,3s', 'T'*27]
        ]

      paths  =  orfer.find_orfs_from_problems(problems, :min_orf_length => 0)
      orfer.orf_sequences_from_trails(paths.trails, 0).should == [
        ['(1s:6),2s,(3s:10)', 'ATGGAAAAAAAAAAAAAAAAAAAATAA'],
        [',(1s:3)', 'TAA'],
        ['1s,2s,3s', 'T'*30],
        ['1s,2s,3s', 'T'*27],
        ['1s,2s,3s', 'T'*27]
        ]

    end
  end

  describe 'sequence2AA' do
    it 'should return corresponding amino acids for an orf sequence' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      orfer.sequence2AA('GCTGCCGCAGCG').should == 'AAAA'
      orfer.sequence2AA('CGTCGCCGACGGAGAAGG').should == 'RRRRRR'
      orfer.sequence2AA('AATAAC').should == 'NN'
      orfer.sequence2AA('GATGAC').should == 'DD'
      orfer.sequence2AA('TGTTGC').should == 'CC'
      orfer.sequence2AA('CAACAG').should == 'QQ'
      orfer.sequence2AA('GAAGAG').should == 'EE'
      orfer.sequence2AA('GGTGGCGGAGGG').should == 'GGGG'
      orfer.sequence2AA('CATCAC').should == 'HH'
      orfer.sequence2AA('ATTATCATA').should == 'III'
      orfer.sequence2AA('TTATTGCTTCTCCTACTG').should == 'LLLLLL'
      orfer.sequence2AA('AAAAAG').should == 'KK'
      orfer.sequence2AA('ATG').should == 'M'
      orfer.sequence2AA('TTTTTC').should == 'FF'
      orfer.sequence2AA('CCTCCCCCACCG').should == 'PPPP'
      orfer.sequence2AA('TCTTCCTCATCGAGTAGC').should == 'SSSSSS'
      orfer.sequence2AA('ACTACCACAACG').should == 'TTTT'
      orfer.sequence2AA('TGG').should == 'W'
      orfer.sequence2AA('TATTAC').should == 'YY'
      orfer.sequence2AA('GTTGTCGTAGTG').should == 'VVVV'
      lambda { orfer.sequence2AA('TAA') }.should raise_error
      lambda { orfer.sequence2AA('TGA') }.should raise_error
      lambda { orfer.sequence2AA('TAG') }.should raise_error
      lambda { orfer.sequence2AA('ABCXYZ') }.should raise_error
    end
  end
end

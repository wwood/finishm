require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "SingleCoherentPathsBetweenNodes" do
  it 'should validate_last_node_of_path_by_recoherence simple' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2,3]
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    finder.validate_last_node_of_path_by_recoherence(paths[0], 15).should == true
  end

  it 'should not validate_last_node_of_path_by_recoherence due to kmer decoherences' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2],
      [2,3],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    finder.validate_last_node_of_path_by_recoherence(paths[0], 15).should == false
  end

  it 'should validate_last_node_of_path_by_recoherence due to kmer decoherences, but too short kmer recoherence' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2],
      [2,3],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    finder.validate_last_node_of_path_by_recoherence(paths[0], 5).should == true
    finder.validate_last_node_of_path_by_recoherence(paths[0], 10).should == true
    finder.validate_last_node_of_path_by_recoherence(paths[0], 11).should == false
  end

  it 'should find a hello world trail' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2,3],
      ])
    initial_path = GraphTesting.make_onodes(graph, %w(1s))
    terminal_oriented_node = GraphTesting.make_onodes(graph, %w(3s)).trail[0]
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new; recoherence_kmer = 15
    problems = finder.find_all_problems(graph, initial_path, terminal_oriented_node, nil, recoherence_kmer)
    #pp problems
    paths = finder.find_paths_from_problems(problems, recoherence_kmer)
    #pp paths
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3],
    ]
  end

  it 'should not find a uncoherent trail' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2],
      [2,3],
      ])
    initial_path = GraphTesting.make_onodes(graph, %w(1s))
    terminal_oriented_node = GraphTesting.make_onodes(graph, %w(3s)).trail[0]
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new; recoherence_kmer = 15
    problems = finder.find_all_problems(graph, initial_path, terminal_oriented_node, nil, recoherence_kmer)
    #pp problems
    paths = finder.find_paths_from_problems(problems, recoherence_kmer)
    #pp paths
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
    ]
  end

  it 'should find one path when only one is coherent 1' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],
    ], 1, 3)
    GraphTesting.add_noded_reads(graph,[
      [1,2,3],
      [1,5],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 15)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3],
    ]
  end

  it 'should find one path when only one is coherent and two nodes at the end' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],

      [3,4],
    ], 1, 4)
    GraphTesting.add_noded_reads(graph,[
      [1,2,3,4],
      [1,5],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 15)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
    ]
  end

  it 'should find both paths when the recoherence_kmer is not long enough' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],
    ], 1, 3)
    GraphTesting.add_noded_reads(graph,[
      [1,2,3],
      [1,5],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 9)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3],
    ]
  end

  it 'should find paths in a complexish road' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [3,4],
      [4,5],

      [1,6],
      [6,3],
      [3,7],
      [7,5],
    ], 1, 5)
    GraphTesting.add_noded_reads(graph,[
      [1,2,3,4,5],
      [1,6],
      [6,3],
      [3,7],
      [7,5],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 15)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4,5],
    ]
  end

  describe 'sub_kmer_sequence_overlap?' do
    # In the LastGraph file, a read is not associated with a node unless in has an
    # entire kmer in it. However, we want more than this for the purposes of recohering
    # across three nodes, otherwise correct paths fail validation.
    graph = Bio::Velvet::Graph.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','LastGraph'))
    sequences = Bio::Velvet::Sequence.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','Sequences'))
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

    it 'should validate a simple path' do
      finder.sub_kmer_sequence_overlap?(
        GraphTesting.make_onodes(graph, '86s,51s,68e'),
        sequences,
        ).should == true
    end

    it 'should not validate a bad first node' do
      finder.sub_kmer_sequence_overlap?(
        GraphTesting.make_onodes(graph, '42s,51s,68e'),
        sequences,
        ).should == false
    end

    it 'should not validate a bad final node' do
      finder.sub_kmer_sequence_overlap?(
        GraphTesting.make_onodes(graph, '86s,51s,14e'),
        sequences,
        ).should == false
    end

#     it 'should validate a longish path' do
#       finder.validate_last_node_of_path_by_recoherence(
#         GraphTesting.make_onodes(graph, '1s,2s,3s,99e,68s,51e,86e,58e,93s'),
#         51,
#         sequences,
#         ).should == true
#     end
  end
end

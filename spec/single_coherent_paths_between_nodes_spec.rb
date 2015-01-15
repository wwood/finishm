require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
require 'bio-velvet'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "SingleCoherentPathsBetweenNodes" do
  it 'should validate_last_node_of_path_by_recoherence simple' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2,3]
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    finder.validate_last_node_of_path_by_recoherence(paths[0], 15, Bio::Velvet::Sequences.new).should == true
  end

  it 'should correctly validate_last_node_of_path_by_recoherence due to kmer decoherences' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2],
      [2,3],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    finder.validate_last_node_of_path_by_recoherence(paths[0], 15, Bio::Velvet::Sequences.new).should == true #kmer must be > 10+7-1+2=18 for read decoherence to click in
    finder.validate_last_node_of_path_by_recoherence(paths[0], 19, Bio::Velvet::Sequences.new).should == false
    finder.validate_last_node_of_path_by_recoherence(paths[0], 18, Bio::Velvet::Sequences.new).should == false
    finder.validate_last_node_of_path_by_recoherence(paths[0], 17, Bio::Velvet::Sequences.new).should == true
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
    problems = finder.find_all_problems(graph, initial_path, terminal_oriented_node, nil, recoherence_kmer, Bio::Velvet::Sequences.new)
    #pp problems
    paths = finder.find_paths_from_problems(problems, recoherence_kmer)
    #pp paths
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3],
    ]
  end

  it 'should not find an incoherent trail' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2],
      [2,3],
      ])
    initial_path = GraphTesting.make_onodes(graph, %w(1s))
    terminal_oriented_node = GraphTesting.make_onodes(graph, %w(3s)).trail[0]
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new; recoherence_kmer = 20
    problems = finder.find_all_problems(graph, initial_path, terminal_oriented_node, nil, recoherence_kmer, Bio::Velvet::Sequences.new)
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
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
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
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
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
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 9, Bio::Velvet::Sequences.new)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3],
      [1,5,3],
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
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4,5],
    ]
  end

  it 'should use recoherence to get around circuits where possible' do
    fail
  end

  it 'should find paths in a graph with a circuit' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],
      [5,5], #circuit
      [3,4]
    ],1 ,4)
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, nil, nil)
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,5,3,4],
    ]
  end

  it 'should find paths with cycles in a graph with a circuit' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],
      [5,5],
      [3,4]
    ], 1, 4)
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, nil, nil, { max_cycles: 1 })
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,5,3,4],
      [1,5,5,3,4]
    ]
  end

  it 'should find paths with cycles in a graph with a multi-entry/exit circuit' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2], #enter circuit 2-3-5-2
      [2,3],
      [3,5],
      [5,2],
      [1,5], #enter circuit 2-3-5-2
      [2,4], #exit circuit
      [3,4] #exit circuit
    ], 1, 4)
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, nil, nil, { max_cycles: 1 })
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,2,3,5,2,3,4],
      [1,2,3,5,2,4],
      [1,2,4],
      [1,5,2,3,4],
      [1,5,2,3,5,2,3,4],
      [1,5,2,3,5,2,4],
      [1,5,2,4],
    ]
  end

  it 'should find paths with cycles in a graph with a figure 8 circuit' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [3,2],
      [3,5],
      [5,3],
      [5,4]
    ], 1, 4)
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, nil, nil, { max_cycles: 1 })
    paths.circular_paths_detected.should == true
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,2,3,5,4],
      [1,2,3,2,3,5,3,2,3,5,4],
      [1,2,3,2,3,5,3,5,4],
      [1,2,3,5,4],
      [1,2,3,5,3,2,3,5,4],
      [1,2,3,5,3,2,3,5,3,5,4],
      [1,2,3,5,3,5,4],
    ]
  end

  it 'should not be subject duplication that may comes from leash length problems' do
    graph, initial_path, terminal_onode = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,5],
      [5,3],
      [3,4]
    ], 1, 4)
    GraphTesting.add_noded_reads(graph,[
      [1,2,3,4],
      [1,5,3,4],
      ])
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,5,3,4],
    ]

    # Make each path shorter in turn
    graph.nodes[5].ends_of_kmers_of_node = 'A'*8
    graph.nodes[5].ends_of_kmers_of_twin_node = 'T'*8
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,5,3,4],
    ]

    # and the other
    graph.nodes[2].ends_of_kmers_of_node = 'A'*7
    graph.nodes[2].ends_of_kmers_of_twin_node = 'T'*7
    paths = finder.find_all_connections_between_two_nodes(graph, initial_path, terminal_onode, nil, 18, Bio::Velvet::Sequences.new)
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
      [1,2,3,4],
      [1,5,3,4],
    ]
  end

  it 'should respect max_explore_nodes' do
    graph, paths = GraphTesting.emit_otrails([
      [1,2,3],
    ])
    GraphTesting.add_noded_reads(graph,[
      [1,2,3],
      ])
    initial_path = GraphTesting.make_onodes(graph, %w(1s))
    terminal_oriented_node = GraphTesting.make_onodes(graph, %w(3s)).trail[0]
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new; recoherence_kmer = 15
    problems = finder.find_all_problems(graph, initial_path, terminal_oriented_node, nil, recoherence_kmer, Bio::Velvet::Sequences.new,
      :max_explore_nodes => 1
      )
    #pp problems
    paths = finder.find_paths_from_problems(problems, recoherence_kmer)
    #pp paths
    paths.circular_paths_detected.should == false
    GraphTesting.sorted_paths(paths.trails).should == [
    ]
  end

  describe 'validate paths by recoherence' do
    graph = Bio::Velvet::Graph.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','LastGraph'))
    sequences = Bio::Velvet::Sequences.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','Sequences'))
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

    it 'should not call badness when there is insufficient read length to validate' do
      fail
    end

    it 'should not call badness when there is no reads but the nodes are arranged in serial order' do
      fail
    end

    it 'should validate a longish path' do
      #TODO not really a thought out test just yet
      alltrail = GraphTesting.make_onodes(graph, '1s,2s,3s,99e,68s,51e,86e,58e,93s')
      otrail = Bio::Velvet::Graph::OrientedNodeTrail.new
      [51,66,77,85].each do |recoherence_kmer|
        (0...alltrail.trail.length).each do |fin|
          otrail.trail = alltrail.trail[0..fin]

          finder.validate_last_node_of_path_by_recoherence(
            otrail,
            recoherence_kmer,
            sequences
            ).should == true
        end
      end
    end
  end

  describe 'sub_kmer_sequence_overlap?' do
    # In the LastGraph file, a read is not associated with a node unless in has an
    # entire kmer in it. However, we want more than this for the purposes of recohering
    # across three nodes, otherwise correct paths fail validation.
    graph = Bio::Velvet::Graph.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','LastGraph'))
    sequences = Bio::Velvet::Sequences.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','Sequences'))
    finder = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

    fwd_reads_graph = Bio::Velvet::Graph.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','LastGraph'))
    fwd_reads_graph.nodes.each do |node|
      node.short_reads.reject!{|r| r.direction == false}
    end
    rev_reads_graph = Bio::Velvet::Graph.parse_from_file(File.join(TEST_DATA_DIR,'gapfilling','5','velvet51_3.5','LastGraph'))
    rev_reads_graph.nodes.each do |node|
      node.short_reads.reject!{|r| r.direction == true}
    end

    it 'should validate a simple path' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(graph, '94s,95s,89s'), sequences).should == true
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '94s,95s,89s'), sequences).should == true
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '94s,95s,89s'), sequences).should == true
    end

    it 'should not validate a bad first node' do
      finder.sub_kmer_sequence_overlap?(
        GraphTesting.make_onodes(graph, '42s,51s,68e'),
        sequences
        ).should == false
    end

    it 'should validate a path where the middle node is reverse' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '4s,46e,93e'), sequences).should == true
    end

    it 'should validate a path where the first node is reverse' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '27e,99e,14s'), sequences).should == true
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '27e,99e,14s'), sequences).should == true
    end

    it 'should validate a path where the last node is reverse' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '86s,51s,68e'), sequences).should == true
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '86s,51s,68e'), sequences).should == true
    end

    it 'should not validate a bad path where the middle node is reverse' do
      # Correct sequence is 4s,46e,93e
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '4s,46e,38s'), sequences).should == false
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '4s,46e,38s'), sequences).should == false
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '20s,46e,93e'), sequences).should == false
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '20s,46e,93e'), sequences).should == false
    end

    it 'should not validate a bad path where the first node is reverse' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '27e,99e,68s'), sequences).should == false
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '27e,99e,68s'), sequences).should == false
    end

    it 'should not validate a bad path where the last node is reverse' do
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(fwd_reads_graph, '86s,51s,14e'), sequences).should == false
      finder.sub_kmer_sequence_overlap?(GraphTesting.make_onodes(rev_reads_graph, '86s,51s,14e'), sequences).should == false
    end

    it 'should validate properly when there is >3 nodes being validated' do
      fail "need to make another assembly for this"
    end

    it 'should validate when multiple nodes at the end do not have that read listed' do
      #e.g. when there is a 1bp node 2nd last, and the last node's length is at least 2bp shorter than the kmer
      # not needed I don't think?
    end
  end

  describe 'CycleFromStartCounter' do
    counter = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder::CycleFromStartCounter.new

    it 'should check for a minimum number of cycles in the start of a sequence of ids' do
      counter.starts_with_minimum_repeats_of_size([1,2,3], 1, 1).should == true
      counter.starts_with_minimum_repeats_of_size([1,2,3], 1, 2).should == false
      counter.starts_with_minimum_repeats_of_size([1,2,1,2,1,2,3], 2, 2).should == true
      counter.starts_with_minimum_repeats_of_size([1,2,1,2,1,2,3], 2, 3).should == true
      counter.starts_with_minimum_repeats_of_size([1,2,1,2,1,2,3], 2, 4).should == false
      counter.starts_with_minimum_repeats_of_size([1,1,2,1,1,2,1,1,2,1,2], 1, 2).should == true
      counter.starts_with_minimum_repeats_of_size([1,1,2,1,1,2,1,1,2,1,2], 1, 3).should == false
      counter.starts_with_minimum_repeats_of_size([1,1,2,1,1,2,1,1,2,1,2], 3, 2).should == true
      counter.starts_with_minimum_repeats_of_size([1,1,2,1,1,2,1,1,2,1,2], 3, 4).should == false
      counter.starts_with_minimum_repeats_of_size([], 4, 1).should == false
      counter.starts_with_minimum_repeats([1,2,3], 1).should == [true, 1]
      counter.starts_with_minimum_repeats([1,1,2,1,1,2,1,1,2,1,2], 3).should == [true, 3]
      counter.starts_with_minimum_repeats([1,1,2,1,1,2,1,1,2,1,2], 4).should == [false, nil]
    end
  end
end

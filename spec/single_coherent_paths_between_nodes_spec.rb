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

#   it 'should calculate a very straightforward trail' do
#     graph, initial_path, terminal_node = GraphTesting.emit_ss([
#       [1,2],
#       [2,3],
#     ], 1, 3)
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes(graph, initial_path, terminal_node, nil)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3],
#     ]
#   end

#   it 'should calculate a trail with just one bubble' do
#     graph, initial_path, terminal_node = GraphTesting.emit_ss([
#       [1,2],
#       [1,3],
#       [2,4],
#       [3,4],
#       [4,5],
#     ], 1, 5)
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes(graph, initial_path, terminal_node, nil)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,4,5],
#       [1,3,4,5],
#     ]
#   end

#   it 'should deal with cycles' do
#     graph, initial_path, terminal_node = GraphTesting.emit_ss([
#       [1,2],
#       [2,3],
#       [3,2],
#       [3,4],
#     ], 1, 4)
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes(graph, initial_path, terminal_node, nil)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4],
#     ]
#   end

#   it 'should find paths not both ending at terminal node' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,3],
#       [3,4],
#       [1,5],
#       [5,3]
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[4]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4],
#       [1,5,3,4],
#     ]
#   end

#   it 'should find through consecutive loops ending at terminal' do
#     # 1 2/3 4 5/6 7
#     graph = GraphTesting.emit([
#       [1,2],
#       [1,3],
#       [2,4],
#       [3,4],
#       [4,5],
#       [4,6],
#       [5,7],
#       [6,7],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[7]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,4,5,7],
#       [1,2,4,6,7],
#       [1,3,4,5,7],
#       [1,3,4,6,7],
#     ]
#   end

#     it 'should find through consecutive loops not ending at terminal' do
#     # 1 2/3 4 5/6 7 8
#     graph = GraphTesting.emit([
#       [1,2],
#       [1,3],
#       [2,4],
#       [3,4],
#       [4,5],
#       [4,6],
#       [5,7],
#       [6,7],
#       [7,8]
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[8]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,4,5,7,8],
#       [1,2,4,6,7,8],
#       [1,3,4,5,7,8],
#       [1,3,4,6,7,8],
#     ]
#   end

#   it 'should find loop off loop 1' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,3],
#       [2,4],
#       [3,8],
#       [4,5],
#       [4,6],
#       [5,7],
#       [6,7],
#       [7,8],
#       [8,9],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[9]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,8,9],
#       [1,2,4,5,7,8,9],
#       [1,2,4,6,7,8,9],
#     ].sort
#   end

#   it 'should find loop off loop with three way' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [1,3],
#       [1,4],
#       [2,5],
#       [3,5],
#       [4,8],
#       [5,6],
#       [5,7],
#       [6,9],
#       [7,9],
#       [8,9],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[9]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,5,6,9],
#       [1,2,5,7,9],
#       [1,3,5,6,9],
#       [1,3,5,7,9],
#       [1,4,8,9],
#     ].sort
#   end

#   it 'should not fail when there is one path beyond the leash (1st), and another not (2nd)' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,3],
#       [1,3],
#       [3,4],
#     ])
#     graph.hash_length = 87
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[4]
#     graph.nodes[1].ends_of_kmers_of_node = 'A'*10
#     graph.nodes[2].ends_of_kmers_of_node = 'A'*100
#     graph.nodes[3].ends_of_kmers_of_node = 'A'*10
#     graph.nodes[4].ends_of_kmers_of_node = 'A'*10
#     (1..4).each do |node_id|
#       graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#     end
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 200, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4],
#       [1,3,4],
#     ].sort

#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 100, true)
#     (1..4).each do |node_id|
#       graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#     end
#     GraphTesting.sorted_paths(paths).should == [
#       [1,3,4],
#     ].sort

#   end


#   it 'probably fails with the nasty leash bug, which is hard to fix' do
#     graph = GraphTesting.emit([
#       [1,3],
#       [1,2],
#       [2,3],
#       [3,4],
#       [4,5],
#     ])
#     graph.hash_length = 87
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[5]
#     graph.nodes[1].ends_of_kmers_of_node = 'A'*10
#     graph.nodes[2].ends_of_kmers_of_node = 'A'*70
#     graph.nodes[3].ends_of_kmers_of_node = 'A'*15
#     graph.nodes[4].ends_of_kmers_of_node = 'A'*10
#     graph.nodes[5].ends_of_kmers_of_node = 'A'*10
#     (1..4).each do |node_id|
#       graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#     end
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 200, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4,5],
#       [1,3,4,5],
#     ].sort

#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 80, true)
#     (1..4).each do |node_id|
#       graph.nodes[node_id].ends_of_kmers_of_twin_node = graph.nodes[node_id].ends_of_kmers_of_node
#     end
#     GraphTesting.sorted_paths(paths).should == [
#       [1,3,4,5],
#     ].sort

#   end

#   it 'should not fail in this special case I realised might trip up the algorithm' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,3],
#       [3,4],
#       [4,5],

#       [1,6],
#       [6,3],
#       [6,5],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[5]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 99999, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4,5],
#       #[1,2,3,4,6,5],
#       [1,6,3,4,5],
#       [1,6,5],
#     ].sort
#   end

#   it 'should not fail in another special case I realised might trip up the algorithm' do
#     #NOTE: to fix this, one must first fix the above graph problem. Argh.
#     # The problem is that a simple joining of the golden path 1,2,3,4,5 and the
#     # golden fragment 1,6,3 yields a circular path
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,3],
#       [3,4],
#       [4,5],

#       [1,6],
#       [6,3],
#       [6,5],

#       [4,6],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[5]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, 99999, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3,4,5],
#       [1,2,3,4,6,5],
#       [1,6,3,4,5],
#       [1,6,5],
#     ].sort
#   end

#   it 'should give the same answer as a more straightfoward (naive?) repeated depth first search style' do
#     raise "need to do some simulation work here to write the test"
#   end

#   it 'should not get confused by a 1 node cycle' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,2],
#       [2,3],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[3]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3],
#     ].sort
#   end

#   it 'should not get confused by a 2 node cycle' do
#     graph = GraphTesting.emit([
#       [1,2],
#       [2,4],
#       [4,2],
#       [2,3],
#     ])
#     initial_node = graph.nodes[1]
#     terminal_node = graph.nodes[3]
#     cartographer = Bio::AssemblyGraphAlgorithms::PathsBetweenNodesFinder.new
#     paths = cartographer.find_all_connections_between_two_nodes_adapter(graph, initial_node, terminal_node, nil, true)
#     GraphTesting.sorted_paths(paths).should == [
#       [1,2,3],
#     ].sort
#   end
end

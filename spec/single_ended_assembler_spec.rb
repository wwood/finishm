require 'ds'
require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "SingleEndedAssembler" do
  describe 'is_short_tip?' do
    it 'should clip short tips correctly' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new(graph)
      assembler.assembly_options[:max_tip_length] = 35
      assembler.is_short_tip?(initial_path[0]).should == [true, [[2, :start_is_first], [3,:start_is_first]]]
      assembler.is_short_tip?(initial_path[0])[0].should == true
      assembler.assembly_options[:max_tip_length] = 25
      assembler.is_short_tip?(initial_path[0])[0].should == false
    end

    it 'should clip short tips in a harder situation' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [3,4],
        [1,4]
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 35
      assembler.is_short_tip?(initial_path[0])[0].should == false
      assembler.assembly_options[:max_tip_length] = 45
      assembler.is_short_tip?(initial_path[0])[0].should == true
    end

    it 'should clip when there is only 1 node' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 25
      assembler.is_short_tip?(initial_path[0])[0].should == true
      assembler.assembly_options[:max_tip_length] = 15
      assembler.is_short_tip?(initial_path[0])[0].should == false
    end
  end




  describe 'assemble_from' do
    it 'should hello world' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s,3s'
    end

    it 'should ignore short tips' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,4],
        [4,5],
        [5,6],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 15
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s,4s,5s,6s'
    end

    it 'should stop at a proper fork' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,4],
        [4,5],
        [5,6],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 5
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s'
    end

    it 'should use single read decoherence' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,4]
        ], 1)
      GraphTesting.add_noded_reads(graph,[
        [1,2,4],
        [2,3]
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 5
      assembler.assembly_options[:recoherence_kmer] = 22
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s,4s'
    end

    it 'should fail when single read decoherence fails' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,4]
        ], 1)
      GraphTesting.add_noded_reads(graph,[
        [1,2],
        [2,3],
        [2,4],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 5
      assembler.assembly_options[:recoherence_kmer] = 22
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s'
    end

    it 'should fail when still forked when single read decoherence fails' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,4]
        ], 1)
      GraphTesting.add_noded_reads(graph,[
        [1,2,3],
        [1,2,4],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 5
      assembler.assembly_options[:recoherence_kmer] = 22
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s'
    end

    it 'should recognise circles' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [3,1],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s,3s,1s'
    end


    it 'should respect the leash length' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [3,1],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:leash_length] = 15
      observed, visits = assembler.assemble_from(initial_path)
      observed.to_shorthand.should == '1s,2s'
    end
  end




  describe 'assemble' do
    it 'should hello world' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 2
      assembler.assembly_options[:min_contig_size] = 0
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      paths.collect{|path| path.to_shorthand}.should == [
        '1s,2s,3s'
        ]
    end

    it 'should assemble two disconnected components' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],

        [4,5],
        [5,6],
        [6,7],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 2
      assembler.assembly_options[:min_contig_size] = 0
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      paths.collect{|path| path.to_shorthand}.should == [
        '1s,2s,3s',
        '4s,5s,6s,7s',
        ]
    end

    it 'should yield correctly' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],

        [4,5],
        [5,6],
        [6,7],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 2
      assembler.assembly_options[:min_contig_size] = 0
      expecteds = DS::Queue.new
      expecteds.enqueue '1s,2s,3s'
      expecteds.enqueue '4s,5s,6s,7s'
      assembler.assemble do |path|
        path.to_shorthand.should == expecteds.dequeue
      end
    end

    it 'should respect the minimum starting node coverage option' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],

        [4,5],
        [5,6],
        [6,7],
        ])
      graph.nodes.each do |node|
        node.coverages = [10]
      end
      graph.nodes[1].coverages[0] = 100
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 2
      assembler.assembly_options[:min_contig_size] = 0
      assembler.assembly_options[:min_coverage_of_start_nodes] = 2

      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      paths.collect{|path| path.to_shorthand}.should == [
        '1s,2s,3s',
        ]
    end

    it 'should not start assembling from short tips' do
      graph = GraphTesting.emit([
        [1,4], #this is the short tip
        [2,3],
        [3,4],
        [4,5],
        [5,6],
        [6,7],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 15
      assembler.assembly_options[:min_contig_size] = 0
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      GraphTesting.sorted_fwd_shorthand_paths(paths).should == [
        '2s,3s,4s,5s,6s,7s',
        ]
    end

    it 'should not start assembling from complex short tips' do
      graph = GraphTesting.emit([
        [1,10], #this is the short tip
        [10,11],
        [10,12],
        [11,6],
        [12,6],

        [2,3],
        [3,4],
        [4,5],
        [5,6],
        [6,7],
        [7,8],
        [8,9],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = 35
      assembler.assembly_options[:min_contig_size] = 40
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      GraphTesting.sorted_fwd_shorthand_paths(paths).should == [
        '2s,3s,4s,5s,6s,7s,8s,9s',
        ]
    end

    it 'should not throw out large single node paths' do
      graph = GraphTesting.emit([
        [1,2],
        ])
      graph.delete_nodes_if{|node| node.node_id == 2}
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = -1
      assembler.assembly_options[:min_contig_size] = 0
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      paths.collect{|path| path.to_shorthand}.should == [
        '1s'
        ]
    end

    it 'should not choke when the starting node is itself forked' do
      graph = GraphTesting.emit([
        [1,2],
        [3,1],
        [3,4],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = -1
      assembler.assembly_options[:min_contig_size] = 5
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      GraphTesting.sorted_fwd_shorthand_paths(paths).should == [
        "2e,1e,3e",
        '3s,4s',
        ]
    end

    it 'should not choke when two nodes along the starting trail are legitimate forks' do
      # a nasty one. Starts at node 1, which then leads to starting at node 3.
      # both node 3 and node 2 have in-arcs
      graph = GraphTesting.emit([
        [3,2],
        [2,1],
        [7,2],
        [6,3],
        [3,8],
        [3,4],
        [3,5],
        ])
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new graph
      assembler.assembly_options[:max_tip_length] = -1
      assembler.assembly_options[:min_contig_size] = 5
      paths = assembler.assemble
      paths.kind_of?(Array).should == true
      GraphTesting.sorted_fwd_shorthand_paths(paths).should == [
        '3e,2e,1e',
        '3e,2e,7e',
        '3s,4s',
        '3s,5s',
        '6s,3s',
        '3s,8s',
        ]
    end

    it 'should be able to deal with recoherence when getting out of a short tip' do
      raise
    end
  end
end

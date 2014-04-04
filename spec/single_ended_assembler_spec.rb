require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "SingleEndedAssembler" do
  describe 'is_short_tip?' do
    it 'should clip short tips correctly' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      assembler.is_short_tip?(initial_path[0], graph, 35).should == true
      assembler.is_short_tip?(initial_path[0], graph, 25).should == false
    end

    it 'should clip short tips in a harder situation' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [3,4],
        [1,4]
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      assembler.is_short_tip?(initial_path[0], graph, 35).should == false
      assembler.is_short_tip?(initial_path[0], graph, 45).should == true
    end

    it 'should clip when there is only 1 node' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      assembler.is_short_tip?(initial_path[0], graph, 25).should == true
      assembler.is_short_tip?(initial_path[0], graph, 15).should == false
    end
  end

  describe 'assemble_from' do
    it 'should hello world' do
      graph, initial_path = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        ], 1)
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      observed = assembler.assemble_from(initial_path, graph)
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
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      observed = assembler.assemble_from(initial_path, graph, :max_tip_length => 15)
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
      assembler = Bio::AssemblyGraphAlgorithms::SingleEndedAssembler.new
      observed = assembler.assemble_from(initial_path, graph, :max_tip_length => 5)
      observed.to_shorthand.should == '1s,2s'
    end
  end
end

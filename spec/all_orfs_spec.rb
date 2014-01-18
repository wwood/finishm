require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')


describe "AllOrfs" do
  describe 'orfs_from_start_stop_indices' do
    it 'should work when there are no orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([],[],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[],[],[]]
      res.final_start_positions.should == [nil]*3
    end

    it 'should work for one orf' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0],[6],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[[0,6]],[],[]]
      res.final_start_positions.should == [nil]*3
    end

    it 'should work for one orf in 2 frames' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0,2],[6,11],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[[0,6]],[],[[2,11]]]
      res.final_start_positions.should == [nil]*3
    end

    it 'should work for one orf and a leftover' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([0,9],[6],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[[0,6]],[],[]]
      res.final_start_positions.should == [9,nil,nil]
    end

    it 'should work with a start orf' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([7],[],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[],[],[]]
      res.final_start_positions.should == [nil,7,nil]
    end

    it 'should work with 3 orfs' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[11,17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[],[],[[8,11],[14,17],[20,44]]]
      res.final_start_positions.should == [nil,nil,nil]
    end

    it 'should work with an internal start codon' do
      orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
      res = orfer.orfs_from_start_stop_indices([8,14,20],[17,44],0)
      res.kind_of?(Bio::AssemblyGraphAlgorithms::AllOrfsFinder::ORFsResult).should == true
      res.start_stop_pairs.should == [[],[],[[8,17],[20,44]]]
      res.final_start_positions.should == [nil,nil,nil]
    end
  end
end

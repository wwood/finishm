require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

class DummyTrail
  attr_accessor :sequence

  def initialize(seq)
    @sequence = seq
  end

  def self.trails(seqs)
    seqs.collect{|s| DummyTrail.new(s)}
  end
end

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "kmer coverage path filter" do
  it 'should rule things out when no appropriate kmers are present' do
    kmers = Bio::KmerMultipleAbundanceHash.new
    kmers['AAAA'] = [0, 0]
    trail = DummyTrail.new('ATATATTTA')
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [1]*2)
    paths.should == []
  end

  it 'should not rule things out when appropriate kmers are present' do
    kmers = Bio::KmerMultipleAbundanceHash.new
    kmers['ATATATTT'] = [1, 10]
    kmers['TATATTTA'] = [3, 1]
    trail = DummyTrail.new('ATATATTTA')
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [1]*2)
    paths.collect{|t| t.sequence}.should == %w(ATATATTTA)
  end

  it 'should not accept different thresholds for different timepoints' do
    kmers = Bio::KmerMultipleAbundanceHash.new
    kmers['ATATATTT'] = [1, 10]
    kmers['TATATTTA'] = [1, 3]
    trail = DummyTrail.new('ATATATTTA')
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [2]*2)
    paths.should == []
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [0,2])
    paths.collect{|t| t.sequence}.should == %w(ATATATTTA)
  end

  it 'should rule some things out in a figure 8 style graph' do
    kmers = Bio::KmerMultipleAbundanceHash.new
    kmers['ATA'] = [1, 10]
    kmers['AGA'] = [1, 10]
    kmers['TAG'] = [1, 10]
    kmers['GAT'] = [1, 10]
    kmers['ATC'] = [1, 10]
    kmers['AGC'] = [1, 10]
    trails = DummyTrail.trails %w(ATATC ATAGC AGATC AGAGC)
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter(trails, kmers, [1]*2)
    paths.collect{|t| t.sequence}.should == %w(ATAGC AGATC)
  end

  it 'should respect exclusion of filtering at the ends' do
    kmers = Bio::KmerMultipleAbundanceHash.new
    kmers['ATATATTT'] = [1, 10]
    kmers['TATATTTA'] = [3, 1]
    trail = DummyTrail.new('GATATATTTAC')
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [1]*2)
    paths.collect{|t| t.sequence}.should == %w()
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [1]*2, :exclude_ending_length => 1)
    paths.collect{|t| t.sequence}.should == %w(GATATATTTAC)
    trail = DummyTrail.new('GATATATTTAGC')
    paths = Bio::AssemblyGraphAlgorithms::KmerCoverageBasedPathFilter.new.filter([trail], kmers, [1]*2, :exclude_ending_length => 1)
    paths.collect{|t| t.sequence}.should == %w()
  end
end

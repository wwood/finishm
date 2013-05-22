require 'rspec'
require 'pp'
require 'systemu'
require 'spec_helper'

require 'kmer_abundance_pattern'


describe 'kmer_abundance_pattern' do
  it 'should parse from human' do
    pat = KmerAbundancePattern.new
    pat.length.should eq(0)

    pat.parse_from_human '0111'
    pat.should eq([false, true, true, true])

    pat.parse_from_human '101'
    pat.should eq([true, false, true])
  end

  it 'should same_as?' do
    pat = KmerAbundancePattern.new
    exp = KmerAbundancePattern.new
    pat.parse_from_human '101'

    exp.parse_from_human '101'
    pat.same_as?(exp).should eq(true)

    exp.parse_from_human '111'
    pat.same_as?(exp).should eq(false)

    exp.parse_from_human '110'
    pat.same_as?(exp).should eq(false)
  end

  it 'should consistent_with?' do
    pat = KmerAbundancePattern.new
    exp = KmerAbundancePattern.new
    pat.parse_from_human '101'

    exp.parse_from_human '101'
    pat.consistent_with?(exp).should eq(true)

    exp.parse_from_human '111'
    pat.consistent_with?(exp).should eq(true)

    exp.parse_from_human '110'
    pat.consistent_with?(exp).should eq(false)

    exp.parse_from_human '11-'
    pat.consistent_with?(exp).should eq(true)

    exp.parse_from_human '1-0'
    pat.consistent_with?(exp).should eq(false)
  end
end

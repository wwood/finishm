require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')


describe "ReadInput" do
  it 'should break contigs easy' do
    breaker = Bio::FinishM::ScaffoldBreaker.new

    Tempfile.open('a') do |tmp|
      tmp.puts '>a'
      tmp.puts 'AAAAA'
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes.length.should == 1
      s = brokes[0]
      s.kind_of?(Bio::FinishM::ScaffoldBreaker::Scaffold).should == true
      s.name.should == 'a'
      s.contigs.length.should == 1
      c = s.contigs[0]
      c.scaffold.should == s
      c.scaffold_position_start.should == 1
      c.scaffold_position_end.should == 5
      c.length.should == 5
    end
  end


  it 'should break contigs 2' do
    breaker = Bio::FinishM::ScaffoldBreaker.new

    Tempfile.open('a') do |tmp|
      tmp.puts '>ab'
      tmp.puts 'AAAAANNNGGG'
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes.length.should == 1
      s = brokes[0]
      s.kind_of?(Bio::FinishM::ScaffoldBreaker::Scaffold).should == true
      s.name.should == 'ab'
      s.contigs.length.should == 2

      c = s.contigs[0]
      c.scaffold.should == s
      c.scaffold_position_start.should == 1
      c.scaffold_position_end.should == 5
      c.length.should == 5

      c = s.contigs[1]
      c.scaffold.should == s
      c.scaffold_position_start.should == 9
      c.scaffold_position_end.should == 11
      c.length.should == 3
    end
  end

  it 'should break contigs 2' do
    breaker = Bio::FinishM::ScaffoldBreaker.new

    Tempfile.open('a') do |tmp|
      tmp.puts '>ab'
      tmp.puts 'AAAAANNNGGG'
      tmp.puts '>bab'
      tmp.puts 'NNAAAAANNNGGG'
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes.length.should == 2
      s = brokes[1]
      s.kind_of?(Bio::FinishM::ScaffoldBreaker::Scaffold).should == true
      s.name.should == 'bab'
      s.contigs.length.should == 2

      c = s.contigs[0]
      c.scaffold.should == s
      c.scaffold_position_start.should == 3
      c.scaffold_position_end.should == 7
      c.length.should == 5

      c = s.contigs[1]
      c.scaffold.should == s
      c.scaffold_position_start.should == 11
      c.scaffold_position_end.should == 13
      c.length.should == 3
    end
  end

  it 'should give back good gap objects' do
     breaker = Bio::FinishM::ScaffoldBreaker.new

    Tempfile.open('a') do |tmp|
      tmp.puts '>ab'
      tmp.puts 'AAAAANNNGGG'
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes.length.should == 1
      gaps = brokes[0].gaps
      gaps.length.should == 1
      gaps[0].start.should == 6
      gaps[0].stop.should == 8
      gaps[0].number.should == 0
    end
  end

  it 'should give back good 3 gap objects' do
     breaker = Bio::FinishM::ScaffoldBreaker.new

    Tempfile.open('a') do |tmp|
      tmp.puts '>ab'
      tmp.puts 'AAAAANNNGGGNNTTNNAA'
      #         1234567890123456789
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes.length.should == 1
      gaps = brokes[0].gaps
      gaps.length.should == 3
      gaps[0].start.should == 6
      gaps[0].stop.should == 8
      gaps[0].number.should == 0
      gaps[1].start.should == 12
      gaps[1].stop.should == 13
      gaps[1].number.should == 1
      gaps[2].start.should == 16
      gaps[2].stop.should == 17
      gaps[2].number.should == 2
    end
  end

  it 'should give back the correct sequence' do
    breaker = Bio::FinishM::ScaffoldBreaker.new
    Tempfile.open('a') do |tmp|
      tmp.puts '>ab'
      seq = 'AAAAANNNGGGNNTTNNAA'
      tmp.puts seq
      #         1234567890123456789
      tmp.close

      brokes = breaker.break_scaffolds(tmp.path)
      brokes[0].sequence.should == seq
    end
  end
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio-commandeer'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "FinishM assemble" do
  it 'should assemble something easy' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    trails = nil
    Dir.chdir(data_path) do
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --quiet --output-contigs /dev/stdout --assemble-from '1e' --fasta 2seqs.sammy.fa"#, :log => log
    end
    splits = trails.split("\n")
    splits.length.should == 2

    splits[0].should == ">1e"
    splits[1][0..60].should == 'TATCCGGTCCCCCTAGAATGTTGATTCCTCCGTCTTCTGATTTCCGTTGGCGGTTCGTATC'
    splits[1].length.should == 613
  end

  it 'should output a pathspec' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    trails = nil
    Dir.chdir(data_path) do
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --quiet --output-pathspec --output-contigs /dev/stdout --assemble-from '1e' --fasta 2seqs.sammy.fa"#, :log => log
    end
    splits = trails.split("\n")
    splits.length.should == 2

    splits[0].should == ">1e 1e"
    splits[1][0..60].should == 'TATCCGGTCCCCCTAGAATGTTGATTCCTCCGTCTTCTGATTTCCGTTGGCGGTTCGTATC'
    splits[1].length.should == 613
  end

  it 'should work when assembling the entire graph' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    trails = nil
    Dir.chdir(data_path) do
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --no-progressbar --quiet --output-pathspec --output-contigs /dev/stdout --fasta 2seqs.sammy.fa"#, :log => log
    end
    splits = trails.split("\n")
    splits.length.should == 6

    splits[0].should == ">contig1 1s"
    splits[1][0..60].should == 'AGGGCAGATTCCCACGCGTTACGCACCCGTGCGCCACTAGACCCGAAGGTCTCGTTCGACT'
    splits[1].length.should == 613
  end

  it 'should not crash when recoherence kmer is given' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    trails = nil
    Dir.chdir(data_path) do
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --no-progressbar --assembly-kmer 7 --recoherence-kmer 51 --quiet --output-pathspec --output-contigs /dev/stdout --fasta 2seqs.sammy.fa"#, :log => log
    end
    # Not sure what is supposed to come out of here with such a short assembly-kmer, but just so long as it doesn't crash.
  end
end

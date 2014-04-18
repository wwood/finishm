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
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --no-progressbar --quiet --output-pathspec --min-contig-length 0 --output-contigs /dev/stdout --fasta 2seqs.sammy.fa"#, :log => log
    end
    splits = trails.split("\n")
    splits.length.should == 24

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

  it 'should output contig stats' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    stats = nil
    Dir.chdir(data_path) do
      stats = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} assemble --no-progressbar --min-contig-length 0 --quiet --output-pathspec --output-contigs /dev/null --output-contig-stats /dev/stdout --fasta 2seqs.sammy.fa"#, :log => log
    end
    #        +contig1 21.15370018975332
    #        +contig2 62.07488986784141
    #        +contig3 22.654135338345863
    #        +contig4 34.90963855421687
    # ...
    stats.match(/name\tcoverage\ncontig1\t21.15370018975332\ncontig2\t62.07488986784141\ncontig3\t22.654135338345863\n/).nil?.should == false
  end
end

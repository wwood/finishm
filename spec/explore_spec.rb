require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio-commandeer'

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "FinishM explore" do
  it 'should explore just a single but rather complex trail' do
    data_path = File.join(TEST_DATA_DIR,'explore','1')
    trails = nil
    Dir.chdir(data_path) do
      trails = Bio::Commandeer.run "#{FINISHM_SCRIPT_PATH} explore --quiet --contigs a.fa --interesting-ends random:end --output-explored-paths - --fasta 2seqs.sammy.fa"
    end
    splits = trails.split("\n").to_a.sort
    splits.length.should == 128
    splits[0].match(/>random:end Dead end \/ coverage nodes:(\d),[\d,]+/)[1] == '1'
    splits[127][0...100].should == 'CCTATCCGGTCCCCCTAGAATGTTGATTCCTCCGTCTTCTGATTTCCGTTGGCGGTTCGTATCGGCTCCCGTACAACCGCCGGAATTGAAGTGTACCTTG'
  end
end

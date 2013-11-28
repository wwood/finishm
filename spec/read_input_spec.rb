require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')


describe "ReadInput" do
  it 'should output correct fasta velvet singles' do
    input = Bio::FinishM::ReadInput.new
    input.fasta_singles = ['my.fasta']
    input.velvet_read_arguments.should == ' -fasta -short my.fasta'
    input.fasta_singles = ['my.fasta','another.fasta']
    input.velvet_read_arguments.should == ' -fasta -short my.fasta -short2 another.fasta'
  end
end

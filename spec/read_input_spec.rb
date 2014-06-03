require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'tempfile'

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')


describe "ReadInput" do
  it 'should output correct fasta velvet singles' do
    input = Bio::FinishM::ReadInput.new
    input.fasta_singles = ['my.fasta']
    input.velvet_read_arguments.should == ' -fasta -short my.fasta'
    input.fasta_singles = ['my.fasta','another.fasta']
    input.velvet_read_arguments.should == ' -fasta -short my.fasta another.fasta'
  end

#   it 'should handle interleaved pairs' do
#     input = Bio::FinishM::ReadInput.new
#     input.interleaved_fastq = ['my.fastq']
#     input.velvet_read_arguments.should == ' -fastq -shortPaired my.fastq'
#     input.interleaved_fastq_gz = ['my.fastq','fq2']
#     input.velvet_read_arguments.should == ' -fastq -shortPaired my.fastq -fastq.gz -shortPaired my.fastq fq2'
#   end
end

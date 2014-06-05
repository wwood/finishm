require 'tempfile'
require 'bio'
require 'bio-commandeer'
require 'spec_helper'

describe 'finishm wander' do
  path_to_script = File.join(File.dirname(__FILE__),'..','bin','finishm wander')

  it 'should scripting test ok with a 1 node thing' do
    random = nil
    Bio::FlatFile.foreach("#{TEST_DATA_DIR }/wander/1/random1.fa") do |s|
      random = s.seq
      break
    end

    Tempfile.open('testing') do |t|
    #File.open('/tmp/contigs','w') do |t|
      t.puts '>first300'
      t.puts random[0...300]
      t.puts '>last400'
      t.puts random[600..-1]
      t.close
      command = "#{path_to_script} --quiet --fasta #{TEST_DATA_DIR}/wander/1/random1.sammy.fa --contigs #{t.path} --output-connections /dev/stdout --overhang 100 --assembly-kmer 51"
      #puts command
      output = Bio::Commandeer.run command
      puts output
      output.split("\n").should == [
        "first300:end\tlast400:start\t500"
        ] #TODO: fix this bug where the distance is incorrectly reported
    end
  end

  it 'should work like scripting 1 node thing except calling it directly' do
    random = nil
    Bio::FlatFile.foreach("#{TEST_DATA_DIR }/wander/1/random1.fa") do |s|
      random = s.seq
      break
    end

    Tempfile.open('testing') do |t|
      #File.open('/tmp/contigs','w') do |t|
      t.puts '>first300'
      t.puts random[0...300]
      t.puts '>last400'
      t.puts random[600..-1]
      t.close

      Tempfile.open('output_connections') do |conns|
        conns.close

        wanderer = Bio::FinishM::Wanderer.new
        wanderer.run(Bio::FinishM::Wanderer::DEFAULT_OPTIONS.merge(
          Bio::FinishM::GraphGenerator::DEFAULT_OPTIONS).merge({
            :contigs_file => t.path,
            :output_connection_file => conns.path,
            :fasta_singles => "#{TEST_DATA_DIR}/wander/1/random1.sammy.fa",
            :contig_end_length => 100,
            }))
      end
    end
  end

  it 'should work with kmer recoherence' do
    random = nil
    Bio::FlatFile.foreach("#{TEST_DATA_DIR }/wander/1/random1.fa") do |s|
      random = s.seq
      break
    end

    Tempfile.open('testing') do |t|
      #File.open('/tmp/contigs','w') do |t|
      t.puts '>first300'
      t.puts random[0...300]
      t.puts '>last400'
      t.puts random[600..-1]
      t.close

      Tempfile.open('output_connections') do |conns|
        conns.close

        wanderer = Bio::FinishM::Wanderer.new
        wanderer.run(Bio::FinishM::Wanderer::DEFAULT_OPTIONS.merge(
          Bio::FinishM::GraphGenerator::DEFAULT_OPTIONS).merge({
            :contigs_file => t.path,
            :output_connection_file => conns.path,
            :fasta_singles => "#{TEST_DATA_DIR}/wander/1/random1.sammy.fa",
            :contig_end_length => 100,
            :recoherence_kmer => 71,
            }))
      end
    end
  end
end

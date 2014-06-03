require 'spec_helper'

describe 'graph_generator' do
  it 'should run using minimal options' do
    probes = [
      'AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC',
      'ACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGCGTAACGCGTATGCAATCTGCCTTGTACTAAGGGATAGCCCAGAGA',
      ]
    read_inputs = Bio::FinishM::ReadInput.new
    read_inputs.fasta_singles_gz = [
      File.join(File.dirname(__FILE__),'data','gapfilling','3','reads.fa.gz')
      ]
    probed_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probes, read_inputs, {
      :assembly_coverage_cutoff => 0,
      :velvet_kmer_size => 31,
      })

    expect(probed_graph).to be_kind_of(Bio::FinishM::ProbedGraph)
  end

  it 'should parse the sequence file or not
    Dir.mktmpdir do |tmpdir|
      probed_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probes, read_inputs, {
        :assembly_coverage_cutoff => 0,
        :velvet_kmer_size => 31,
        })

      expect(probed_graph).to be_kind_of(Bio::FinishM::ProbedGraph)
    end
end

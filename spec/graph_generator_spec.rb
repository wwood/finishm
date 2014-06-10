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
    expect(probed_graph.velvet_sequences).to be_kind_of(Bio::Velvet::Underground::BinarySequenceStore)
    expect(probed_graph.velvet_sequences[1]).to eq('AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC')
    expect(probed_graph.velvet_sequences[2]).to eq('ACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGCGTAACGCGTATGCAATCTGCCTTGTACTAAGGGATAGCCCAGAGA')
  end

  it 'should check probe sequences are present in the assembly' do
    Dir.mktmpdir do |tmpdir|
      probes = [
        'AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC',
        'ACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGCGTAACGCGTATGCAATCTGCCTTGTACTAAGGGATAGCCCAGAGA',
        ]
      read_inputs = Bio::FinishM::ReadInput.new
      read_inputs.fasta_singles_gz = [
        File.join(File.dirname(__FILE__),'data','gapfilling','3','reads.fa.gz')
        ]
      # First assembly run is fine
      probed_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probes, read_inputs, {
        :assembly_coverage_cutoff => 0,
        :velvet_kmer_size => 31,
        :output_assembly_path => tmpdir,
        })
      expect(probed_graph).to be_kind_of(Bio::FinishM::ProbedGraph)
      expect(probed_graph.velvet_sequences[1]).to eq('AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC')

      # Reading in the last assembly with the same probes should work
      probed_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probes, read_inputs, {
        :previous_assembly => tmpdir,
        })
      expect(probed_graph).to be_kind_of(Bio::FinishM::ProbedGraph)
      expect(probed_graph.velvet_sequences[1]).to eq('AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC')

      # Reading assembly expecting a different probe read should fail
      expect {Bio::FinishM::GraphGenerator.new.generate_graph(['AAAAAAAAAAA'], read_inputs, {
        :previous_assembly => tmpdir,
        })}.to raise_exception

      expect {Bio::FinishM::GraphGenerator.new.generate_graph(['AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGC',
        'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        ], read_inputs, {
        :previous_assembly => tmpdir,
        })}.to raise_error
    end
  end

  it 'should find probes by probe names' do
    # >32_1   63      0
    # >32_2   64      0
    # >33_1   65      0
    # >33_2   66      0
    # >34_1   67      0
    # >34_2   68      0
    # >35_1   69      0
    # >35_2   70      0
    names = %w(32_1 35_2)

    Dir.mktmpdir do |tmpdir|
      read_inputs = Bio::FinishM::ReadInput.new
      read_inputs.fasta_singles_gz = [
        File.join(File.dirname(__FILE__),'data','gapfilling','3','reads.fa.gz')
        ]
      # First assembly run is fine
      probed_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_inputs, {
        :assembly_coverage_cutoff => 0,
        :velvet_kmer_size => 31,
        :output_assembly_path => tmpdir,
        :probe_read_names => names,
        })
      probed_graph.probe_node_reads.collect{|read| read.read_id}.should == [63,70]
    end
  end

  describe 'check_probe_sequences' do
    it 'should check_probe_sequences' do
      graph_gen = Bio::FinishM::GraphGenerator.new
      expect(graph_gen.check_probe_sequences(['AAA'],['','AAA'])).to eq(true)
      expect(graph_gen.check_probe_sequences(['AAA'],['','ATA'])).to eq(false)
      expect(graph_gen.check_probe_sequences(['TTTA','AAA'],['','TTTA','AAA'])).to eq(true)
      expect(graph_gen.check_probe_sequences(['TTTA','AAA'],['','TTTA','ATA'])).to eq(false)
    end
  end
end

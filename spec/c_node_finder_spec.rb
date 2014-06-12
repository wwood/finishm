require 'spec_helper'

describe 'c_probe_node_finder' do
  it 'should find set of probe nodes' do
    Dir.mktmpdir do |tmpdir|
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
        #:output_assembly_path => '/tmp/v',
        :output_assembly_path => tmpdir,
        })

      finder = Bio::FinishM::CProbeNodeFinder.new
      read_probing_graph = Bio::Velvet::Underground::Graph.parse_from_file File.join(TEST_DATA_DIR,'c_probe_node_finder','1','LastGraph')
      finder.find_probe_nodes(read_probing_graph, [4,5,6]).should == [1]
      finder.find_probe_nodes(read_probing_graph, [515,135]).should == [3,4]
    end
  end

  it 'should return the same as the Ruby node_finder' do
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

      finder = Bio::FinishM::CProbeNodeFinder.new
      read_probing_graph = Bio::Velvet::Underground::Graph.parse_from_file File.join(TEST_DATA_DIR,'c_probe_node_finder','1','LastGraph')
      probes = finder.find_probes(read_probing_graph, [1,544])
      probes[0][0].node_id.should == 1
      probes[1][0].node_id.should == 2
      probes[0][1].should == true
      probes[1][1].should == false
      probes[0][2].read_id.should == 1
      probes[1][2].read_id.should == 544

      finder.find_probes(read_probing_graph, []).should == []
    end
  end
end

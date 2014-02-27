require 'bio-commandeer'

describe 'finishm gap closer' do
  path_to_script = File.join(File.dirname(__FILE__),'..','bin','finishm sequence')
  it 'should scripting test ok' do
    command = "#{path_to_script} --quiet --already-assembled-velvet-directory spec/data/contig_printer/1/seq.fa.velvet --path 9s,12s,7e"
    expected = 'TTAGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTACGGAATAACAGTTAGAAATGACTGCTAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCC'+"\n"
    observed = Bio::Commandeer.run command
    observed.should == expected
  end


  it 'should should fail when bad input is given bad node id' do
    command = "#{path_to_script} --quiet --already-assembled-velvet-directory spec/data/contig_printer/1/seq.fa.velvet --path 9s,12s,7e,0e"
    expect {
      Bio::Commandeer.run command
      }.to raise_error
  end

  it 'should should fail when bad input is given bad dir' do
    command = "#{path_to_script} --quiet --already-assembled-velvet-directory spec/data/contig_printer/1/seq.fa.velvet --path 9z,12s,7e"
    expect {
      Bio::Commandeer.run command
      }.to raise_error
  end

  it 'should should fail when bad input is given bad connection' do
    command = "#{path_to_script} --quiet --already-assembled-velvet-directory spec/data/contig_printer/1/seq.fa.velvet --path 9s,7e"
    expect {
      Bio::Commandeer.run command
      }.to raise_error
  end
end

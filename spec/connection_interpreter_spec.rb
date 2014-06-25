require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.create_connection(first_second,distance=10)
    doing_first = true
    probes = first_second.collect do |namer|
      if namer.kind_of?(Fixnum)
        probe = Bio::FinishM::ConnectionInterpreter::Probe.new
        probe.sequence_name = namer.to_s
        if doing_first
          probe.side = :end
          doing_first = false
        else
          probe.side = :start
        end
        probe #'return'

      elsif matches = namer.match(/^(.+)([se])$/)
        probe = Bio::FinishM::ConnectionInterpreter::Probe.new
        probe.sequence_name = matches[1]
        probe.side = matches[2] == 's' ? :start : :end
        probe #'return'
      else
        raise namer
      end
    end
    conn = Bio::FinishM::ConnectionInterpreter::Connection.new
    conn.probe1 = probes[0]
    conn.probe2 = probes[1]
    conn.distance = distance
    return conn
  end

  def self.create_connections(array_of_conn_strings)
    conns = array_of_conn_strings.collect{|a| create_connection a}
    seqs = {}
    conns.each_with_index do |conn, i|
      [conn.probe1, conn.probe2].each_with_index do |probe, j|
        seqs[probe.sequence_name] ||= (['A']*10).join
      end
    end
    return conns, seqs
  end
end

describe "ConnectionInterpreter" do
  it 'should find doubly_single_contig_connections hello world' do
    conns, seqs = GraphTesting.create_connections([
      %w(contig1s contig3e)
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).doubly_single_contig_connections.collect{|c| c.to_s}.should == [
        'contig1s/contig3e:10'
        ]
  end

  it 'should find inter-contig connections in a loop' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      [3,1],
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).doubly_single_contig_connections.collect{|c| c.to_s}.sort.should == [
        '1e/2s:10',
        '2e/3s:10',
        '3e/1s:10',
        ].sort
  end

  it 'should not include all connections in inter-contig connections' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      [3,4],
      [3,5],
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).doubly_single_contig_connections.collect{|c| c.to_s}.sort.should == [
        '1e/2s:10',
        '2e/3s:10',
        ].sort
  end

  it 'should be able to handle loops where there is only one contig' do
    conns, seqs = GraphTesting.create_connections([
      [1,1],
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).doubly_single_contig_connections.collect{|c| c.to_s}.sort.should == [
        '1e/1s:10',
        ].sort
  end

  it 'should scaffold hello world' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      ])
    observed = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).scaffolds
    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.name.should == 'scaffold1'
    o.sequence.should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.scaffold_position_start}.should == [0,20]
    o.contigs.collect{|c| c.scaffold_position_end}.should == [9,29]
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold 3 contigs together' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      ])
    observed = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).scaffolds
    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.name.should == 'scaffold1'
    o.sequence.should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.scaffold_position_start}.should == [0,20,40]
    o.contigs.collect{|c| c.scaffold_position_end}.should == [9,29,49]
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold two separate scaffolds and a leftover' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [3,4],
      ])
    seqs['contig99'] = 'ATGC'
    observed = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).scaffolds
    observed.should be_kind_of(Array)
    observed.length.should == 3
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.name.should == 'scaffold1'
    o.sequence.should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.scaffold_position_start}.should == [0,20]
    o.contigs.collect{|c| c.scaffold_position_end}.should == [9,29]
    o = observed[1]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.name.should == 'scaffold2'
    o.sequence.should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.scaffold_position_start}.should == [0,20]
    o.contigs.collect{|c| c.scaffold_position_end}.should == [9,29]
    o = observed[2]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.contigs[0].original_name.should == 'contig99'
    o.sequence.should == 'ATGC'
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold single contig circular scaffolds' do
    conns, seqs = GraphTesting.create_connections([
      [1,1],
      ])
    observed = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).scaffolds
    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.name.should == 'singleton1'
    o.sequence.should == 'AAAAAAAAAA'
    o.contigs.collect{|c| c.scaffold_position_start}.should == [0]
    o.contigs.collect{|c| c.scaffold_position_end}.should == [9]
    o.circular?.should == true
  end

  it 'should scaffold multi-contig circular scaffolds' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      [3,1],
      [9,10],
      ])
    seqs['contig99'] = 'ATGC'
    observed = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).scaffolds
    observed.length.should == 3
    observed.collect{|o| o.circular?}.should == [true, false, false]
    observed[0].contigs.collect{|c| c.original_name}.should == %w(3 1 2)
  end
end

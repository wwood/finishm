require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

class GraphTesting
  def self.create_connection(first_second,distance=10)
    doing_first = true
    probes = first_second.collect do |namer|
      if namer.kind_of?(Fixnum)
        probe = Bio::FinishM::ConnectionInterpreter::Probe.new
        probe.sequence_index = namer
        if doing_first
          probe.side = :end
          doing_first = false
        else
          probe.side = :start
        end
        probe #'return'

      elsif matches = namer.match(/^(.+)([se])$/)
        probe = Bio::FinishM::ConnectionInterpreter::Probe.new
        probe.sequence_index = matches[1].to_i
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
        seqs[probe.sequence_index] ||= (['A']*10).join
      end
    end
    return conns, seqs
  end
end

describe "ConnectionInterpreter" do
  it 'should find doubly_single_contig_connections hello world' do
    conns, seqs = GraphTesting.create_connections([
      %w(1s 3e)
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs
      ).doubly_single_contig_connections.collect{|c| c.to_s}.should == [
        '1s/3e:10'
        ]
  end

  it 'should find inter-contig connections in a loop' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      [3,1],
      ])
    Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
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
      conns, seqs.keys
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
      conns, seqs.keys
      ).doubly_single_contig_connections.collect{|c| c.to_s}.sort.should == [
        '1e/1s:10',
        ].sort
  end

  it 'should scaffold hello world' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      ])
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)
    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.contigs.collect{|c| c.sequence_index}.should == [1,2]
    o.contigs.collect{|c| c.direction}.should == [true, true]
    o.gap_lengths.should === [10]
    o.sequence(seqs).should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold 3 contigs together' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      ])
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)
    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.sequence(seqs).should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.sequence_index}.should == [1,2,3]
    o.contigs.collect{|c| c.direction}.should == [true, true, true]
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold two separate scaffolds and a leftover' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [3,4],
      ])
    seqs[99] = 'ATGC'
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)

    observed.should be_kind_of(Array)
    observed.length.should == 3
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.sequence(seqs).should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.sequence_index}.should == [1,2]
    o.contigs.collect{|c| c.direction}.should == [true, true]
    o = observed[1]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.sequence(seqs).should == 'AAAAAAAAAANNNNNNNNNNAAAAAAAAAA'
    o.contigs.collect{|c| c.sequence_index}.should == [3,4]
    o.contigs.collect{|c| c.direction}.should == [true, true]
    o = observed[2]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.contigs[0].sequence_index.should == 99
    o.sequence(seqs).should == 'ATGC'
    observed.collect{|o| o.circular?}.uniq.should == [false]
  end

  it 'should scaffold single contig circular scaffolds' do
    conns, seqs = GraphTesting.create_connections([
      [1,1],
      ])
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)

    observed.should be_kind_of(Array)
    observed.length.should == 1
    o = observed[0]
    o.should be_kind_of(Bio::FinishM::ConnectionInterpreter::Scaffold)
    o.sequence(seqs).should == 'AAAAAAAAAA'
    o.contigs.collect{|c| c.sequence_index}.should == [1]
    o.contigs.collect{|c| c.direction}.should == [true]
    o.circular?.should == true
  end

  it 'should scaffold multi-contig circular scaffolds' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      [2,3],
      [3,1],
      [9,10],
      ])
    seqs[87] = 'ATGC'
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)

    observed.length.should == 3
    observed.collect{|o| o.circular?}.should == [true, false, false]
    observed[0].contigs.collect{|c| c.sequence_index}.should == [3,1,2]
  end

  it 'should respect the given distance given in the Connection' do
    conns, seqs = GraphTesting.create_connections([
      %w(1s 3e)
      ])
    conns[0].distance = 5
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)
    observed.length.should == 1
    observed[0].gap_lengths.should == [5]
  end

  it 'should be able to handle reverse scaffolding' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      ])
    conns[0].probe1.side = :start #reverse both contigs
    conns[0].probe2.side = :end
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, seqs.keys
      )
    observed = interpreter.scaffolds(interpreter.doubly_single_contig_connections)
    observed.length.should == 1
    observed[0].sequence(seqs).should == 'TTTTTTTTTTNNNNNNNNNNTTTTTTTTTT'
  end

  it 'should report unconnected probes' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      ])
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2]
      )
    interpreter.unconnected_probes.collect{|pro| pro.to_settable}.should == [
      [1, :start],
      [2, :end],
      ]
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2,3]
      )
    interpreter.unconnected_probes.collect{|pro| pro.to_settable}.should == [
      [1, :start],
      [2, :end],
      [3, :start],
      [3, :end]
      ]
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2,4]
      )
    interpreter.unconnected_probes.collect{|pro| pro.to_settable}.should == [
      [1, :start],
      [2, :end],
      [4, :start],
      [4, :end]
      ]
  end

  it 'should report unconnected sequences' do
    conns, seqs = GraphTesting.create_connections([
      [1,2],
      ])
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2]
      )
    interpreter.unconnected_sequences.should == [
      ]
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2,3]
      )
    interpreter.unconnected_sequences.should == [
      3
      ]
    interpreter = Bio::FinishM::ConnectionInterpreter.new(
      conns, [1,2,4,5,6]
      )
    interpreter.unconnected_sequences.should == [
      4,5,6
      ]
  end
end

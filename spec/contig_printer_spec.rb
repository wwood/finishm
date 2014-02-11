require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm'); Bio::Log::CLI.configure('bio-velvet')

class Bio::AssemblyGraphAlgorithms::ContigPrinter::Connection
  def comparable
    [
      reference_path.collect{|onode| onode.to_settable}.flatten,
      variants.collect do |variant|
        [variant.to_settable]
      end.sort
      ]
  end
end

describe "ContigPrinter" do
  describe "two_contigs_and_connection_to_printable_connection" do

    it 'should work with just a straight ref path' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3], []
        )
    end

    it 'should work with just 1 variant' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3],
        [1,4,3],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3], [[1,3,4]]
        )
    end

    it 'should work with 3 variants' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3],
        [1,4,3],
        [1,5,3],
        [1,6,3],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3], [
          [1,3,4],
          [1,3,5],
          [1,3,6],
          ])
    end

    it 'should work with a >1 node variant' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3],
        [1,4,5,6,3],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3], [[1,3,4,5,6]]
        )
    end

    it 'should work with 2 bubbles' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3,4,5],
        [1,6,3,7,5],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3,4,5], [[1,3,6],[3,5,7]]
        )
    end

    it 'should work with 2 overlapping bubbles' do
      graph, paths = GraphTesting.emit_paths([
        [1,2,3,4,5],
        [1,6,3,7,5],
        [1,2,3,7,5],
        ])
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      conn = printer.two_contigs_and_connection_to_printable_connection(paths)
      conn.comparable.should == GraphTesting.emit_printer_connection(graph,
        [1,2,3,4,5], [[1,3,6],[3,5,7]]
        )
    end
  end

  describe 'one_connection_between_two_contigs first' do
    it 'should work with a path with two nodes both in the same direction' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_read_id = 161 #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_read_id = 1045
      acon.start_probe_node = graph.nodes[9]
      acon.end_probe_node = graph.nodes[4]
      acon.start_probe_contig_offset = 0
      acon.end_probe_contig_offset = 0
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '12345'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
        '67890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890'
        )
      observed.should == expected
    end

    it 'should handle reads not starting right at the end of the contig' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_read_id = 161 #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_read_id = 1045
      acon.start_probe_node = graph.nodes[9]
      acon.end_probe_node = graph.nodes[4]
      acon.start_probe_contig_offset = 2
      acon.end_probe_contig_offset = 3
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '123'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
        '90'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890'
        )
      observed.should == expected
    end

    it 'should handle reads not starting at the start of the node' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_read_id = 709 #355_1
      acon.end_probe_read_id = 671 #336_1
      acon.start_probe_node = graph.nodes[9]
      acon.end_probe_node = graph.nodes[4]
      acon.start_probe_contig_offset = 1
      acon.end_probe_contig_offset = 1
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '1234'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.
          gsub(/.{6}$/,'').
          gsub(/^.{11}/,'') +
        '7890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890'
        )
      observed.should == expected
    end



#     it 'should work when the direction of the start read is rev' do
#       graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
#       graph.nodes.length.should == 13
#       acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
#       acon.start_probe_read_id = 60 #30_2 #Found these by using bwa and inspecting the Sequence velvet file
#       acon.end_probe_read_id = 1045
#       acon.start_probe_node = graph.nodes[9]
#       acon.end_probe_node = graph.nodes[4]
#       acon.start_probe_contig_offset = 0
#       acon.end_probe_contig_offset = 0
#       acon.paths = [
#         GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
#         ]
#       expected = '12345'+
#         File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
#         '67890'
#       observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
#         graph,'12345',acon,'67890'
#         )
#       observed.should == expected

#     end
  end
end

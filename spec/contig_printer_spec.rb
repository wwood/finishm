require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

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
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "ContigPrinter" do
  describe 'sequences_to_variants_conservative' do
    it 'should handle a multi-variant' do
      seqs = [
              'ATGAATATGTGCATAGGATT',
              'ATGAATCGATGCAGTTGATT',
              #      ###    ###
              #01234567890123456789
             ]
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      ref, variants = printer.send(:sequences_to_variants_conservative, seqs)

      ref.should == 'ATGAATNNNTGCANNNGATT'
      variants.collect{|v| v.to_shorthand
      }.sort.should == [
                   '7S:ATG',
                   '7S:CGA',
                   '14S:TAG',
                   '14S:GTT',
                  ].sort
    end

    it 'should handle a semi-redundant multi-variant' do
      seqs = [
              'ATGAATATGTGCATAGGATT',
              'ATGAATCGATGCAGTTGATT',
              'ATGAATCGATGCATAGGATT',
              #      ###    ###
              #01234567890123456789
             ]
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      ref, variants = printer.send(:sequences_to_variants_conservative, seqs)

      ref.should == 'ATGAATNNNTGCANNNGATT'
      variants.collect{|v| v.to_shorthand
      }.sort.should == [
                   '7S:ATG',
                   '7S:CGA',
                   '14S:TAG',
                   '14S:GTT',
                  ].sort
    end

    it 'should handle not variants' do
      seqs = [
              'ATGAATATGTGCATAGGATT',
             ]
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      ref, variants = printer.send(:sequences_to_variants_conservative, seqs)

      ref.should == 'ATGAATATGTGCATAGGATT'
      variants.collect{|v| v.to_shorthand
      }.sort.should == [].sort
    end

    it 'should handle gaps' do
      seqs = [
              'ATTCTGAACGTAAGCATTATATGAATATGTGCATAGGATTTATTGGATCAGTGGCACGTA',
              'ATTCTGAACGTAAGCATTATATGAATATGTGCAGTTGATTTATTGGATCAGTGGCACGTA',
              'ATTCTGAACGTAAGCATTATATGAATCGATGAGTT GATTTATTGGATCAGTGGCACGTA',
              #                          ###  - !!!
              #1234567890123456789012345678901234567890123456..........
              #         1         2         3         4
             ]
      printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
      ref, variants = printer.send(:sequences_to_variants_conservative, seqs)

      ref.should == 'ATTCTGAACGTAAGCATTATATGAATNNNTGNNNNNGATTTATTGGATCAGTGGCACGTA'
      variants.collect{|v| v.to_shorthand
      }.sort.should == [
        '27S:ATG',
        '27S:CGA',
        '32S:CATAG',
        '32S:CAGTT',
        '32S:AGTT',
        '36D:1',
        ].sort
    end
  end



  #     it 'should work with 3 variants' do
  #       graph, paths = GraphTesting.emit_paths([
  #         [1,2,3],
  #         [1,4,3],
  #         [1,5,3],
  #         [1,6,3],
  #         ])
  #       printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
  #       conn = printer.two_contigs_and_connection_to_printable_connection(paths)
  #       conn.comparable.should == GraphTesting.emit_printer_connection(graph,
  #         [1,2,3], [
  #           [1,3,4],
  #           [1,3,5],
  #           [1,3,6],
  #           ])
  #     end

  #     it 'should work with a >1 node variant' do
  #       graph, paths = GraphTesting.emit_paths([
  #         [1,2,3],
  #         [1,4,5,6,3],
  #         ])
  #       printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
  #       conn = printer.two_contigs_and_connection_to_printable_connection(paths)
  #       conn.comparable.should == GraphTesting.emit_printer_connection(graph,
  #         [1,2,3], [[1,3,4,5,6]]
  #         )
  #     end

  #     it 'should work with 2 bubbles' do
  #       graph, paths = GraphTesting.emit_paths([
  #         [1,2,3,4,5],
  #         [1,6,3,7,5],
  #         ])
  #       printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
  #       conn = printer.two_contigs_and_connection_to_printable_connection(paths)
  #       conn.comparable.should == GraphTesting.emit_printer_connection(graph,
  #         [1,2,3,4,5], [[1,3,6],[3,5,7]]
  #         )
  #     end

  #     it 'should work with 2 overlapping bubbles' do
  #       graph, paths = GraphTesting.emit_paths([
  #         [1,2,3,4,5],
  #         [1,6,3,7,5],
  #         [1,2,3,7,5],
  #         ])
  #       printer = Bio::AssemblyGraphAlgorithms::ContigPrinter.new
  #       conn = printer.two_contigs_and_connection_to_printable_connection(paths)
  #       conn.comparable.should == GraphTesting.emit_printer_connection(graph,
  #         [1,2,3,4,5], [[1,3,6],[3,5,7]]
  #         )
  #     end
  #   end

  describe 'one_connection_between_two_contigs first' do
    it 'should work with a path with two nodes both in the same direction' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_noded_read = graph.nodes[9].short_reads.select{|nr| nr.read_id == 161}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 0
      acon.end_probe_contig_offset = 0
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '12345'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
        '67890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890',[]
        )
      observed.should == expected
    end

    it 'should handle reads not starting right at the end of the contig' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_noded_read = graph.nodes[9].short_reads.select{|nr| nr.read_id == 161}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 1
      acon.end_probe_contig_offset = 2
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '1234'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
        '890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890',[]
        )
      observed.should == expected
    end

    it 'should reverse reads starting from reverse onodes to start the process' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      #  #<Bio::Velvet::Graph::NodedRead:0x0000000293fd08 @direction=false, @offset_from_start_of_node=2, @read_id=289, @start_coord=0>,
      acon.start_probe_noded_read = graph.nodes[7].short_reads.select{|nr| nr.read_id == 289}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 0
      acon.end_probe_contig_offset = 0
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '12345'+
        'TAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCCCAGGGATGAGCCCGCGTAGGATTAGCTTGTTGGTGAGGTAAAGGCTCACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACATGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACTCCGGAATT' +
        '67890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890',[]
        )
      observed.should == expected
    end

    it 'should reverse reads starting from reverse onodes to start the process' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      #  #<Bio::Velvet::Graph::NodedRead:0x0000000293fd08 @direction=false, @offset_from_start_of_node=2, @read_id=289, @start_coord=0>,
      acon.start_probe_noded_read = graph.nodes[7].short_reads.select{|nr| nr.read_id == 289}[0] #Found these by using bwa and inspecting the Sequence velvet file
      #  #<Bio::Velvet::Graph::NodedRead:0x00000002763ef8 @direction=false, @offset_from_start_of_node=3, @read_id=800, @start_coord=0>,
      acon.end_probe_noded_read = graph.nodes[10].short_reads.select{|nr| nr.read_id == 800}[0]
      acon.start_probe_contig_offset = 0
      acon.end_probe_contig_offset = 0
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(7e 13s 5e 11e 2s 10s))
        ]
      expected = '12345'+
        'TAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCCCAGGGATGAGCCCGCGTAGGATTAGCTTGTTGGTGAGGTAAAGGCTCACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACATGGCCCAGACTCCTACGGGAGGCAGCAG' +
        '67890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890',[]
        )
      observed.should == expected
    end
  end

  describe 'ready_two_contigs_and_connections' do
    it 'should handle a single SNP' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_noded_read = graph.nodes[9].short_reads.select{|nr| nr.read_id == 161}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 2
      acon.end_probe_contig_offset = 3
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e)),#highest coverage
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 1e 2e 10s 4e)),
        ]
      expected =
        'ATGAACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTACGG'+
        'AATAACAGTTAGAAATGACTGCTAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCCCAGGGATGAGCCCGCGTAGGATTAGCTTGTTGGTGAGGTAAANN'+
        'NTNNCNNANNNNNNNNNNNNTNNNNNGNNNNNNNNNNNGNTNAGNNNCNNNGNNNNNGNGANNTGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGC'+
        'GAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCCCCGGCTAACTCCGTG'+
        'CCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACTCCGGAATTCA'
      observed_sequence, observed_variants = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.ready_two_contigs_and_connections(
        graph,'ATGCA',acon,'ATGCA',[]
        )
      observed_sequence.should == expected
      observed_variants.collect{|v|v.to_shorthand}.sort.should == [
        "210S:GGC", "214S:CA", "217S:CA", "220S:GGCGACGATCCT", "233S:AGCTG", "239S:TCTGAGAGGAT", "251S:A", "253S:C", "256S:CCA", "260S:ACT", "264S:GGACT", "270S:A", "273S:CA",
        "210S:TTT", "214S:AC", "217S:TC", "220S:CCAACAAGCTAA", "233S:CCTAC", "239S:CGCTCAGACCA", "251S:C", "253S:A", "256S:GAT", "260S:GTC", "264S:CCTTG", "270S:T", "273S:GC",
        ].sort
    end

    it 'should handle when start_coord is not == 0 and both reads are inwards facing' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_noded_read = graph.nodes[9].short_reads.select{|nr| nr.read_id == 161}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 0
      acon.end_probe_contig_offset = 0

      # introduce badness
      acon.start_probe_noded_read.start_coord = 3
      acon.end_probe_noded_read.start_coord = 4
      reads = {
        161 => 'A'*100,
        1045 => 'C'*100,
        }

      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e))
        ]
      expected = '12345'+'AAA'+
        File.open(File.join TEST_DATA_DIR, 'contig_printer','1','seq2_1to550.fa').readlines[1].strip.gsub(/..$/,'') +
        'CCCC'+'67890'
      observed = Bio::AssemblyGraphAlgorithms::ContigPrinter.new.one_connection_between_two_contigs(
        graph,'12345',acon,'67890', reads
        )
      observed.should == expected
    end

    it 'should handle when start_coord is not == 0 and both reads are outwards facing' do
      raise
    end

    it 'should handle when the example path is not the same length as the reference path' do
      fail
    end
  end

  describe 'AnchoredConnection' do
    it 'should collapse_paths_to_maximal_coverage_path!' do
      graph = Bio::Velvet::Graph.parse_from_file(File.join TEST_DATA_DIR, 'contig_printer','1','seq.fa.velvet','LastGraph')
      graph.nodes.length.should == 13
      acon = Bio::AssemblyGraphAlgorithms::ContigPrinter::AnchoredConnection.new
      acon.start_probe_noded_read = graph.nodes[9].short_reads.select{|nr| nr.read_id == 161}[0] #Found these by using bwa and inspecting the Sequence velvet file
      acon.end_probe_noded_read = graph.nodes[4].short_reads.select{|nr| nr.read_id == 1045}[0]
      acon.start_probe_contig_offset = 2
      acon.end_probe_contig_offset = 3
      acon.paths = [
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 11e 2s 10s 4e)),#highest coverage
        GraphTesting.make_onodes(graph, %w(9s 12s 7e 13s 5e 1e 2e 10s 4e)),
        ]
      acon.collapse_paths_to_maximal_coverage_path!
      acon.paths.collect{|path| path.to_shorthand}.should == [%w(9s 12s 7e 13s 5e 11e 2s 10s 4e).join(',')]
    end
  end
end

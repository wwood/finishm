require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

class GraphTesting
  def self.metapath_to_array(metapath)
    to_return = []
    metapath.collect do |node_or_arr|
      if node_or_arr.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
        paths = []
        rev = node_or_arr.is_reverse
        node_or_arr.each_path do |path|
          path_in_bubble = rev ? path[1..-1] : path[0...-1]
          if path_in_bubble.length == 1
            paths.push path_in_bubble[0].node_id
          else
            paths.push path_in_bubble.collect{|onode| onode.node_id}
          end
        end
        if rev
          to_return.push node_or_arr.converging_oriented_node_settable[0]
          to_return.push paths.sort
        else
          to_return.push paths.sort
          to_return.push node_or_arr.converging_oriented_node_settable[0]
        end
      elsif node_or_arr.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        to_return.push node_or_arr.node_id
      else
        raise "Unknown metapath element: #{node_or_arr.inspect}"
      end
    end
    return to_return
  end

  def self.metapaths_to_arrays(metapaths)
    return metapaths.collect do |m|
      metapath_to_array m
    end
  end
end


describe "BubblyAssembler" do
  describe 'assemble_from' do
    it 'should assemble something easy' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        ], 1, 3)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,2,3]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..3).to_a
    end

    it 'should handle the simplest bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,4],
        [1,3],
        [3,4],
        ], 1, 4)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = -1
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,3],4]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..4).to_a
    end

    it 'should handle a simple bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],
        [3,5],
        ], 1, 5)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = -1
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..5).to_a
    end

    it 'should handle a slightly complex bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [3,4],
        [4,7],
        [2,5],
        [5,4],
        [5,6],
        [6,7],
        [7,8]
        ], 1, 8)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = -1
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,2,[[3,4],[5,4],[5,6]],7,8]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..8).to_a
    end

    it 'should deal with fluff in the middle of a bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],
        [3,5],
        [4,99],
        ], 1, 5)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = 11
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..5).to_a+[99]
    end

    it 'should deal with fluff at the beginning of a bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [2,99],
        [1,4],
        [4,3],
        [3,5],
        ], 1, 5)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = 11
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..5).to_a+[99]
    end

    it 'should deal with fluff at the end of a bubble' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],
        [3,5],
        [5,6],
        [3,99]
        ], 1, 6)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = 11
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5,6]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..6).to_a+[99]
    end

    it 'should be able to handle two bubbles' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],
        [3,5],

        [5,6],
        [5,7],
        [6,8],
        [8,9],
        [7,8],
        [7,10],
        [10,9],

        [9,11],
        [11,12],

        [7,99],
        ], 1, 6)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = 11
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5,[[6,8],[7,8],[7,10]],9,11,12]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..12).to_a + [99]
    end

    it 'should respect the leash length' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],

        [3,5],
        [5,6],
        ], 1, 2)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = -1
      cartographer.assembly_options[:leash_length] = 5
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == [1]
    end

    it 'should be able to start from initial paths that include bubbles' do
      graph, initial_path, terminal = GraphTesting.emit_ss([
        [1,2],
        [2,3],
        [1,4],
        [4,3],
        [3,5],
        [5,6],
        [3,99]
        ], 1, 6)
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph
      cartographer.assembly_options[:max_tip_length] = 11
      metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
      GraphTesting.metapath_to_array(metapath).should == [1,[2,4],3,5,6]
      visited_nodes.to_a.collect{|s| s[0]}.sort.should == (1..6).to_a + [99]

      metapath2, visited_nodes2 = cartographer.assemble_from(metapath[0..1], Set.new((1..4)))
      GraphTesting.metapath_to_array(metapath2).should == [1,[2,4],3,5,6]
      visited_nodes2.to_a.collect{|s| s[0]}.sort.should == (1..6).to_a + [99]
    end
  end

  describe 'assemble' do
    it 'should be able to get out of the middle of a simple path' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
        [4,1],
        ])
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
        :max_tip_length => -1,
        :min_contig_size => 0,
        }
      metapaths = cartographer.assemble
      GraphTesting.metapaths_to_arrays(metapaths).should == [[4,1,2,3]]
    end

    it 'should be able to get out of the middle convergent node when both sides have bubbles' do
      graph = GraphTesting.emit([
        [1,2],
        [2,3],
        [4,1],

        [11,12],
        [11,13],
        [12,4],
        [13,4],

        [3,21],
        [3,22],
        [21,23],
        [22,23],
        ])
      cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
        :max_tip_length => -1,
        :min_contig_size => 0,
        }
      metapaths = cartographer.assemble
      GraphTesting.metapaths_to_arrays(metapaths).should == [
        [11,[12,13],4,1,2,3,[21,22],23]
        ]
    end

    it 'should be able to get out of the middle of a path when it is bubbly' do
      fail
    end

    it 'should be able to assemble several contigs' do
      fail
    end

    it 'should be able to use recoherence to get around an inter-genome repeat' do
      fail
    end
  end
end




describe 'metapath' do
  it 'should act like an array' do
    metapath = Bio::AssemblyGraphAlgorithms::BubblyAssembler::MetaPath.new
    metapath.length.should == 0
    graph, initial_path, term = GraphTesting.emit_ss([[1,2]],1,2)
    onode = initial_path[0]
    metapath << onode
    metapath.length.should == 1
    metapath[0].should == onode
  end

  it 'should be able to output a reference contig' do
    graph, initial_path = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      [3,5],
      ],1,1)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
      :max_tip_length => -1,
      :min_contig_size => 0,
      }
    metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
    metapath.reference_trail.to_shorthand.should == '1s,2s,3s,5s'
  end

  it 'should be able to iterate over variants from the reference contig' do
    fail
  end

  it 'should be able to write out a reference and VCF format file' do
    fail
  end
end


describe 'Bubble' do
  it 'should reference trail compare by node ID' do
    graph, initial_path = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      ],1,1)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
      :max_tip_length => -1,
      :min_contig_size => 0,
      }
    metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
    bubble = metapath[1]
    bubble.reference_trail.should be_kind_of(Bio::Velvet::Graph::OrientedNodeTrail)
    bubble.reference_trail.to_shorthand.should == '2s,3s'
  end

  it 'should reference trail compare by coverage' do
    graph, initial_path = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      ],1,1)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
      :max_tip_length => -1,
      :min_contig_size => 0,
      }
    metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
    bubble = metapath[1]
    graph.nodes[2].coverages = [10]
    graph.nodes[4].coverages = [20]
    bubble.reference_trail.to_shorthand.should == '4s,3s'
  end

  it 'should reference trail slightly more complex bubble' do
    graph, initial_path = GraphTesting.emit_ss([
      [2,3],
      [3,4],
      [4,5],
      [2,6],
      [6,4],
      [6,7],
      [7,5],
      ],2,7)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
      :max_tip_length => -1,
      :min_contig_size => 0,
      }
    metapath, visited_nodes = cartographer.assemble_from(initial_path, nil)
    bubble = metapath[1]
    bubble.reference_trail.to_shorthand.should == '3s,4s,5s'
    graph.nodes[4].coverages = [10]
    graph.nodes[7].coverages = [20]
    bubble.reference_trail.to_shorthand.should == '6s,7s,5s'
  end

  it 'should reference trail when path is reversed' do
    graph, initial_path, terminal = GraphTesting.emit_ss([
      [1,2],
      [2,3],
      [1,4],
      [4,3],
      [3,5],
      ], 1, 5)
    cartographer = Bio::AssemblyGraphAlgorithms::BubblyAssembler.new graph, {
      :max_tip_length => -1,
      :min_contig_size => 0,
      }
    metapath, v = cartographer.assemble_from initial_path, nil
    metapath.to_shorthand.should == '1s,{2s,3s|4s,3s},5s'
    metapath.reverse!
    metapath.to_shorthand.should == "5e,{3e,2e|3e,4e},1e"
    metapath.reference_trail.to_shorthand.should == '5e,3e,2e,1e'
  end
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'
class Util
  def self.revcom(seq)
    Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
  end
end

#Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('finishm'); Bio::Log::CLI.configure('finishm')

describe "Fluffer" do
  it 'should work with a very straightforward single path' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
    ], [1])
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 100)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1,2],
      ],
    ]
  end

  it 'should work with two starting points' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,4]
    ], [1,3])
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 100)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1,2,3],
      ],
      [
        [3,4],
      ]
    ]
  end

  it 'should work with two starting points and a dead end' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,4],
      [2,5],
    ], [1,3])
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 100)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1,2,3],
        [1,2,5],
      ],
      [
        [3,4],
      ]
    ]
  end

  it 'should respect leash lengths' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,4],
      [2,5],
    ], [1,3])
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 5)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1],
      ],
      [
        [3],
      ]
    ]
  end

  it 'should be able to handle missing probes' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,4],
      [2,5],
    ], [1,3])
    finishm_graph.probe_nodes[2] = nil
    finishm_graph.probe_node_directions[2] = nil
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 5)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1],
      ],
      [
        [3],
      ],[
      ]
    ]
  end

  it 'should be able to handle missing probes in the middle' do
    finishm_graph = GraphTesting.finishm_graph([
      [1,2],
      [2,3],
      [3,4],
      [2,5],
    ], [1,2,3])
    finishm_graph.probe_nodes[1] = nil
    finishm_graph.probe_node_directions[1] = nil
    fluffer = Bio::AssemblyGraphAlgorithms::Fluffer.new
    paths = fluffer.fluff(finishm_graph, 5)
    GraphTesting.sorted_array_of_paths(paths).should == [
      [
        [1],
      ],
      [
      ],
      [
        [3],
      ]
    ]
  end
end

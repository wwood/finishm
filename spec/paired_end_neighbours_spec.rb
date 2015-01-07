require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "PairedNeighbours" do
  it 'should get a paired neighbour node' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2]])
    GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]
    neighs.collect{|n| n.distance}.should == [80]
    neighs[0].num_adjoining_reads.should == 1
    p neighs[0]
    p neighs[0].node
    neighs[0].node.kind_of?(Bio::Velvet::Graph::Node).should == true
  end

  it 'should get a direct neighbour node' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2]])
    #GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == []
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [2]
    p neighs[0]
    p neighs[0].node
    neighs[0].node.kind_of?(Bio::Velvet::Graph::Node).should == true
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :end_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [1]
  end

  it 'should get both direct and paired neighbours together' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1,1],[2,1],[3,2]])
    GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[1], :start_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [2,3]
    neighs.collect{|n| n.distance}.should == [0,80]
    neighs[1].num_adjoining_reads.should == 1
  end

  it 'should choose the best direction on a paired neighbour node' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1, 1],[2, 1,3,5],[3, 2,4,6]])
    GraphTesting.make_reads_paired(f)

    f.graph.nodes[3].short_reads.find{|r| r.read_id == 6}.offset_from_start_of_node = 5
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]
    neigh = neighs[0]
    neigh.distance.should == 81.66666666666667
    neigh.num_adjoining_reads.should == 3

    # confuse it so that read 6 is now facing a different direction
    f.graph.nodes[3].short_reads.find{|r| r.read_id == 6}.direction = false #previously true
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]
    neigh = neighs[0]
    neigh.distance.should == 80
    neigh.num_adjoining_reads.should == 2
  end

  it 'should work with min number of num_adjoining_reads' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1, 1],[2, 1,3,5],[3, 2,4,6]])
    GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)

    finder.min_adjoining_reads = 0
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]

    finder.min_adjoining_reads = 4
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == []
  end

  it 'should not accept when min_adjoining_reads is violated after orientation checking' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1, 1],[2, 1,3,5],[3, 2,4,6]])
    GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)

    finder.min_adjoining_reads = 3
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]

    f.graph.nodes[3].short_reads.get_read_by_id(4).direction = false
    finder.min_adjoining_reads = 3
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == []
  end

  it 'should avoid high coverage nodes' do
    f = GraphTesting.finishm_graph([[1,2],[3,4]])
    GraphTesting.add_reads_to_nodes(f, [[1, 1],[2, 1,3,5],[3, 2,4,6]])
    GraphTesting.make_reads_paired(f)
    finder = Bio::FinishM::PairedEndNeighbourFinder.new(f, 100)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(f.graph.nodes[2], :start_is_first)

    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]

    f.graph.nodes[3].coverages = [99*10,99*10]

    finder.max_adjoining_node_coverage = 98
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == []

    finder.max_adjoining_node_coverage = 99
    neighs = finder.neighbours(onode)
    neighs.collect{|n| n.node.node_id}.should == [3]
  end
end

require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "DepthFirstSearch" do
  it 'should traverse a linear array' do
    Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug');
    Bio::Log::LoggerPlus.new 'bio-velvet'; Bio::Log::CLI.configure 'bio-velvet'

    graph = Bio::Velvet::Graph.new
    node1 = Bio::Velvet::Graph::Node.new
    node1.node_id = 1
    node2 = Bio::Velvet::Graph::Node.new
    node2.node_id = 2
    node3 = Bio::Velvet::Graph::Node.new
    node3.node_id = 3
    array = Bio::Velvet::Graph::NodeArray.new
    array[1] = node1
    array[2] = node2
    array[3] = node3
    arc1 = Bio::Velvet::Graph::Arc.new
    arc1.begin_node_id = 1
    arc1.end_node_id = 2
    arc1.begin_node_direction = true
    arc1.end_node_direction = true
    arc2 = Bio::Velvet::Graph::Arc.new
    arc2.begin_node_id = 2
    arc2.end_node_id = 3
    arc2.begin_node_direction = true
    arc2.end_node_direction = true
    arc_array = Bio::Velvet::Graph::ArcArray.new
    arc_array.push arc1
    arc_array.push arc2
    graph.nodes = array
    graph.arcs = arc_array


    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node1
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

    expecteds = [1,2,3]
    founds = []

    graph.depth_first_search(onode) do |onode|
      founds.push onode.last.node.node_id
    end
    founds.should == expecteds
  end

  it 'should traverse a bubble' do
    graph = Bio::Velvet::Graph.new
    node1 = Bio::Velvet::Graph::Node.new
    node1.node_id = 1
    node2 = Bio::Velvet::Graph::Node.new
    node2.node_id = 2
    node3 = Bio::Velvet::Graph::Node.new
    node3.node_id = 3
    node4 = Bio::Velvet::Graph::Node.new
    node4.node_id = 4
    array = Bio::Velvet::Graph::NodeArray.new
    array[1] = node1
    array[2] = node2
    array[3] = node3
    array[4] = node4
    graph.nodes = array

    arc_array = Bio::Velvet::Graph::ArcArray.new
    graph.arcs = arc_array
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 2
    arc.end_node_id = 4
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 3
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 3
    arc.end_node_id = 4
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc

    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node1
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

    expecteds = [1,2,4,3]
    founds = []

    graph.depth_first_search(onode) do |onode|
      founds.push onode.last.node.node_id
    end
    founds.should == expecteds
  end

  it 'should deal with cycles' do
    graph = Bio::Velvet::Graph.new
    node1 = Bio::Velvet::Graph::Node.new
    node1.node_id = 1
    node2 = Bio::Velvet::Graph::Node.new
    node2.node_id = 2
    array = Bio::Velvet::Graph::NodeArray.new
    array[1] = node1
    array[2] = node2
    graph.nodes = array

    arc_array = Bio::Velvet::Graph::ArcArray.new
    graph.arcs = arc_array
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 1
    arc.end_node_id = 2
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_id = 2
    arc.end_node_id = 2
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc_array.push arc

    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node1
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

    expecteds = [1,2]
    founds = []

    graph.depth_first_search(onode) do |onode|
      founds.push onode.last.node.node_id
    end
    founds.should == expecteds
  end

  it 'should stop when asked' do
    Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug');
    Bio::Log::LoggerPlus.new 'bio-velvet'; Bio::Log::CLI.configure 'bio-velvet'

    graph = Bio::Velvet::Graph.new
    node1 = Bio::Velvet::Graph::Node.new
    node1.node_id = 1
    node2 = Bio::Velvet::Graph::Node.new
    node2.node_id = 2
    node3 = Bio::Velvet::Graph::Node.new
    node3.node_id = 3
    array = Bio::Velvet::Graph::NodeArray.new
    array[1] = node1
    array[2] = node2
    array[3] = node3
    arc1 = Bio::Velvet::Graph::Arc.new
    arc1.begin_node_id = 1
    arc1.end_node_id = 2
    arc1.begin_node_direction = true
    arc1.end_node_direction = true
    arc2 = Bio::Velvet::Graph::Arc.new
    arc2.begin_node_id = 2
    arc2.end_node_id = 3
    arc2.begin_node_direction = true
    arc2.end_node_direction = true
    arc_array = Bio::Velvet::Graph::ArcArray.new
    arc_array.push arc1
    arc_array.push arc2
    graph.nodes = array
    graph.arcs = arc_array


    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node1
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST

    expecteds = [1]
    founds = []

    graph.depth_first_search(onode) do |onode|
      founds.push onode.last.node.node_id
      false
    end
    founds.should == expecteds
  end
end

require 'yaml'

class Bio::FinishM::BadFormatWriter
  def initialize
    @to_yamlify = []
  end

  def add_metapath(name, metapath)
    to_write = []
    metapath.each do |onode_or_bubble|
      if onode_or_bubble.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        next_to_write = {}
        next_to_write['type'] = 'regular'
        next_to_write['node'] = onode_or_bubble.to_shorthand
        next_to_write['coverage'] = onode_or_bubble.node.coverage
        to_write << next_to_write
      else
        # bubble
        paths = []
        onode_or_bubble.each_path do |path|
          next_to_write = {}
          next_to_write['nodes'] = path.to_shorthand
          next_to_write['coverage'] = path.coverage
          paths << next_to_write
        end
        to_write << {
          'type' => 'bubble',
          'paths' => paths
          }
      end
    end

    @to_yamlify << {
      'contig_name' => name,
      'graph' => to_write.to_yaml
      }
  end

  def yaml
    @to_yamlify.to_yaml
  end

  def write(output_io)
    output_io.print yaml
  end
end
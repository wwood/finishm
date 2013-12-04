require 'set'

module Bio::AssemblyGraphAlgorithms
  class CoverageBasedGraphFilter
    include Bio::FinishM::Logging

    # Remove all nodes from the graph that do not have sufficient coverage
    # (i.e. possibly are sequencing error artefacts)
    #
    # Options:
    # :whitelisted_sequences: provide an enumerable of sequence IDs, don't remove any nodes that have reads tracked to any of these IDs
    #
    # Returns nodes_removed, arcs_removed (as objects, in particular order)
    def remove_low_coverage_nodes(graph, threshold, options = {})
      graph.delete_nodes_if do |node|
        deleting = false
        if node.coverage and (node.coverage < threshold)
          deleting = true
        end

        if deleting and options[:whitelisted_sequences]
          options[:whitelisted_sequences].each do |seq_id|
            if node.short_reads.collect{|r| r.read_id}.include?(seq_id)
              deleting = false
              log.debug "Preserving low coverage but whitelisted node: #{node.node_id}" if log.debug?
            end
          end
        end
        deleting
      end
    end
  end

  class ConnectivityBasedGraphFilter
    include Bio::FinishM::Logging

    # Remove parts of the graph that are unconnected to any whitelisted nodes
    def remove_unconnected_nodes(graph, whitelisted_nodes)
      # Copy the whitelist
      all_whitelisted_nodes = Set.new whitelisted_nodes

      # Depth-first search of all the connected parts looking for nodes to keep
      whitelisted_nodes.each do |originally_whitelisted_node|
        [:start_is_first, :end_is_first].each do |direction|
          onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
          onode.node = originally_whitelisted_node
          onode.first_side = direction
          log.debug "Testing for connectivity from #{onode.node.node_id}/#{onode.first_side}" if log.debug?
          graph.depth_first_search(onode) do |path|
            node = path.last
            log.debug "Whitelisting node #{node.node.node_id}/#{node.first_side}" if log.debug?
            all_whitelisted_nodes << node.node
            true
          end
        end
      end

      # Delete all nodes that aren't in the
      graph.delete_nodes_if do |node|
        !all_whitelisted_nodes.include?(node)
      end
    end
  end
end

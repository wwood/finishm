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
    #
    # options:
    # :leash_length: don't explore more than this length away from each of the whitelisted_nodes. Defualt nil, no bounds
    def remove_unconnected_nodes(graph, whitelisted_nodes, options={})
      # Copy the whitelist
      all_whitelisted_nodes = Set.new whitelisted_nodes

      dij = Bio::AssemblyGraphAlgorithms::Dijkstra.new
      dij_options = {:ignore_directions => true}
      dij_options[:leash_length] = options[:leash_length]

      # Depth-first search of all the connected parts looking for nodes to keep
      whitelisted_nodes.each do |originally_whitelisted_node|
        onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
        onode.node = originally_whitelisted_node
        onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST #irrelevant which is first because :ignore_directions => true
        log.debug "Testing for connectivity from #{onode.node.node_id}" if log.debug?

        min_distances = dij.min_distances(graph, onode, dij_options)
        min_distances.each do |key, distance|
          all_whitelisted_nodes << graph.nodes[key[0]]
        end
      end

      # Delete all nodes that aren't in the
      graph.delete_nodes_if do |node|
        !all_whitelisted_nodes.include?(node)
      end
    end
  end
end

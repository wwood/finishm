module Bio::AssemblyGraphAlgorithms
  class CoverageBasedGraphFilter
    # Remove all nodes from the graph that do not have sufficient coverage
    # (i.e. possibly are sequencing error artefacts)
    #
    # Returns nodes_removed, arcs_removed (as objects, in particular order)
    def remove_low_coverage_nodes(graph, threshold)
      graph.delete_nodes_if do |node|
        node.coverage < threshold
      end
    end
  end
end

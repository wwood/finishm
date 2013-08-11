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
        deleting = (node.coverage < threshold)

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
end

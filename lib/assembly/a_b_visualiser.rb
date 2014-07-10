require 'graphviz'
require 'set'

class Bio::Velvet::Graph::Node
  def includes_kmers?(list_of_kmers)
    list_of_kmers.each do |kmer|
      return true if ends_of_kmers_of_node.include?(kmer) or ends_of_kmers_of_twin_node.include?(kmer)
    end
    return false
  end
end

module Bio
  module Assembly
    class ABVisualiser
      include Bio::FinishM::Logging

      # Visualise a (velvet) graph, as a graphviz object
      #
      # Possible options:
      # :start_kmers: list of kmers to denote the start node(s)
      # :end_kmers: list of kmers to denote the end node(s)
      # :start_node_id: ID of node to mark as a start
      # :end_node_id:ID of node to mark as a end
      # :start_node_ids: array of node IDs to mark as a start
      # :end_node_ids:array of node IDs to mark as a end
      # :coverage_cutoff: ignore nodes with less coverage than this cutoff
      # :digraph: output as a digraph (default true, else output undirected graph)
      # :nodes: an Enumerable of nodes to be visualised.
      def graphviz(graph, options={})
        opts = {}
        opts[:type] = :digraph unless options[:digraph] == false
        opts[:overlap] = :scale
        graphviz = GraphViz.new(:G, opts)

        nodes_to_explore = options[:nodes]
        nodes_to_explore ||= graph.nodes

        # Add all the nodes
        blacklisted_node_ids = Set.new
        log.debug "Converting nodes to GraphViz format"
        nodes_to_explore.each do |node|
          cov = node.coverage
          if options[:coverage_cutoff] and cov < options[:coverage_cutoff] and !cov.nil?
            blacklisted_node_ids.add node.node_id
          else
            cov_string = cov.nil? ? '' : cov.round
            mods = {
              :label => "n#{node.node_id}_length#{node.ends_of_kmers_of_node.length}_coverage#{cov_string}",
            }
            includes_start = false
            includes_end = false
            if options[:start_kmers]
              includes_start = node.includes_kmers?(options[:start_kmers])
            end
            if options[:end_kmers]
              includes_end = node.includes_kmers?(options[:end_kmers])
            end
            if options[:start_node_id]
              includes_start = true if node.node_id == options[:start_node_id]
            end
            if options[:end_node_id]
              includes_end = true if node.node_id == options[:end_node_id]
            end
            if options[:start_node_ids]
              includes_start = true if options[:start_node_ids].include? node.node_id
            end
            if options[:end_node_ids]
              includes_end = true if options[:end_node_ids].include? node.node_id
            end

            if includes_start and includes_end
              log.warn "Start and end kmers detected in the same node!"
            elsif includes_start
              mods[:color] = "red"
            elsif includes_end
              mods[:color] = "green"
            end

            graphviz.add_nodes node.node_id.to_s, mods
          end
        end

        # Add all the edges
        arcs_of_interest = graph.arcs
        if options[:nodes]
          arcs_of_interest = Set.new
          nodes_to_explore.each do |node|
            graph.arcs.get_arcs_by_node_id(node.node_id).each do |arc|
              arcs_of_interest << arc
            end
          end
        end
        log.info "Converting #{arcs_of_interest.length} arcs to GraphViz format"
        arcs_of_interest.each do |arc|
          # Add unless the node has been blacklisted
          unless blacklisted_node_ids.include? arc.begin_node_id or
            blacklisted_node_ids.include? arc.end_node_id

            # Direction of the arrows, to denote connection to beginning of node (connects to start = in-arrow-head to node on output graph)
            if arc.connects_end_to_beginning?(arc.begin_node_id, arc.end_node_id)
              graphviz.add_edges arc.begin_node_id.to_s, arc.end_node_id.to_s
            elsif  arc.connects_end_to_end?(arc.begin_node_id, arc.end_node_id)
              graphviz.add_edges arc.begin_node_id.to_s, arc.end_node_id.to_s, {:arrowhead => "none"}
            elsif  arc.connects_beginning_to_beginning?(arc.begin_node_id, arc.end_node_id)
              graphviz.add_edges arc.begin_node_id.to_s, arc.end_node_id.to_s, {:color => "blue"}
            elsif  arc.connects_beginning_to_end?(arc.begin_node_id, arc.end_node_id)
              graphviz.add_edges arc.end_node_id.to_s, arc.begin_node_id.to_s
            end
          end
        end

        return graphviz
      end


      class SimplifiedGraph
        def self.create_from_velvet_graph(graph)
          nodes_incorporated = 0

          # While there is more of the graph to incorporate
          while nodes_incorporated < graph.nodes.length
            raise "not implemented"
          end
        end

        # A class representing a linear string of nodes without forks
        class Path
        end
      end
    end
  end
end

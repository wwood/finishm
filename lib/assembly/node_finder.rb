module Bio
  module AssemblyGraphAlgorithms
    class NodeFinder
      include Bio::FinishM::Logging

      # Find the node whre the read with the given sequence_id (sequence ID from within velvet)
      # resides. Assumes that read tracking was turned on during velvetg command.
      # Returns [node_id, direction] where direction is true if the read is forward-facing
      # relative to the node, or nil if no node could be found.
      # When multiple nodes contain the read, return the node that is closest to the beginning
      # of the read.
      def find_unique_node_with_sequence_id(graph, sequence_id)
        nodes_with_read = graph.nodes.select do |node|
          node.short_reads.select{|r| r.read_id == sequence_id}.length > 0
        end
        log.debug "Found #{nodes_with_read.length} nodes with the anchor read in it: #{nodes_with_read.collect{|n| n.node_id}.sort.join(',')}"
        return nil if nodes_with_read.empty?

        # TODO: There is a slight bit of imperfection here - multiple nodes can be the minimum
        # Hopefully won't be a bit problem
        best_node = nodes_with_read.min do |n1, n2|
          r1 = n1.short_reads.find{|r| r.read_id == sequence_id}
          r2 = n2.short_reads.find{|r| r.read_id == sequence_id}
          r1.offset_from_start_of_node <=> r2.offset_from_start_of_node
        end
        best_noded_read = best_node.short_reads.find{|r| r.read_id == sequence_id}
        return best_node, best_noded_read.direction
      end

      def find_unique_nodes_with_sequence_ids(graph, sequence_ids)
        # Create data structure
        endings = {}
        sequence_ids.each do |seq_id|
          endings[seq_id] = []
        end

        # Fill data structure with candidate nodes
        graph.nodes.each do |node|
          node.short_reads.each do |read|
            if endings[read.read_id]
              endings[read.read_id].push node
            end
          end
        end

        # Pick the best node from each of the candidate nodes for each sequence_id
        endings.collect do |sequence_id, nodes|
          best_node = nodes.min do |n1, n2|
            r1 = n1.short_reads.find{|r| r.read_id == sequence_id}
            r2 = n2.short_reads.find{|r| r.read_id == sequence_id}
            r1.offset_from_start_of_node <=> r2.offset_from_start_of_node
          end
          if best_node
            best_noded_read = best_node.short_reads.find{|r| r.read_id == sequence_id}
            [best_node, best_noded_read.direction]
          else
            []
          end
        end

        return endings
      end

      def find_unique_node_with_kmers(velvet_graph, kmers)
        # TODO: search in a more sane way, algorithmically
        # TODO: only choose kmers that are unique to the assembly
        found_node = nil
        found_direction = nil
        kmers.each_with_index do |fwd, i|
          rev = Bio::Sequence::NA.new(fwd).reverse_complement.to_s.upcase
          #log.debug "Testing kmer #{i} #{fwd} / #{rev}"
          current_haul = []
          velvet_graph.nodes.each do |node|
            if log.debug?
              str = "Searching for #{fwd} and #{rev} in node #{node.node_id}"
              #str += node.sequence
            end
            if node.sequence? and (node.sequence.include?(fwd) or node.sequence.include?(rev))
              current_haul.push node
            end
          end
          if !current_haul.empty?
            if current_haul.length == 1
              node = current_haul[0]
              log.debug "Found a possibly suitable node in the graph, of length #{node.length}" if log.debug?
              found_node = node
              if found_node.sequence.include?(fwd)
                if found_node.sequence.include?(rev)
                  log.debug "Found a kmer that is included in the forward and reverse directions of the same node. Unlucky, ignoring this kmer" if log.debug?
                else
                  log.info "Found a suitable kmer fwd: #{fwd}"
                  found_direction = true
                end
              else
                log.info "Found a suitable kmer rev: #{rev}"
                found_direction = false
              end
              break unless found_direction.nil?
            else
              log.debug "kmer #{fwd}/#{rev} was found in multiple nodes, so not using it as a starting/ending node" if log.debug?
            end
          end
        end
        return found_node, found_direction
      end

      def find_nodes_with_kmers(velvet_graph, kmers)
        # TODO: search in a more sane way, algorithmically
        # TODO: only choose kmers that are unique to the assembly
        found_nodes = []
        kmers.each_with_index do |fwd, i|
          rev = Bio::Sequence::NA.new(fwd).reverse_complement.to_s.upcase
          #log.debug "Testing kmer #{i} #{fwd} / #{rev}"
          current_haul = []
          velvet_graph.nodes.each do |node|
            if log.debug?
              str = "Searching for #{fwd} and #{rev} in node #{node.node_id}"
              #str += node.sequence
            end
            if node.sequence? and (node.sequence.include?(fwd) or node.sequence.include?(rev))
              current_haul.push node
            end
          end
          if !current_haul.empty?
            if current_haul.length == 1
              node = current_haul[0]
              log.debug "Found a possibly suitable node in the graph, of length #{node.length}" if log.debug?
              found_node = node
              if found_node.sequence.include?(fwd)
                if found_node.sequence.include?(rev)
                  log.debug "Found a kmer that is included in the forward and reverse directions of the same node. Unlucky, ignoring this kmer" if log.debug?
                  found_nodes = [found_node]
                else
                  log.debug "Found a suitable kmer fwd: #{fwd}"
                  found_nodes = [found_node]
                end
              else
                log.debug "Found a suitable kmer rev: #{rev}"
                found_nodes = [found_node]
              end
              break
            else
              log.debug "kmer #{fwd}/#{rev} was found in multiple nodes" if log.debug?
              found_nodes = current_haul
            end
          end
        end
        return found_nodes
      end
    end
  end
end


module Bio
  module AssemblyGraphAlgorithms
    class NodeFinder
      def log
        Bio::Log::LoggerPlus['finishm']
      end

      def find_unique_node_with_kmers(velvet_graph, kmers)
        # TODO: search in a more sane way, algorithmically
        # TODO: only choose kmers that are unique to the assembly
        found_node = nil
        found_direction = nil
        kmers.each do |fwd|
          rev = Bio::Sequence::NA.new(fwd).reverse_complement.to_s.upcase
          current_haul = []
          velvet_graph.nodes.each do |node|
            if log.debug?
              str = "Searching for #{fwd} and #{rev} in node #{node.node_id}/"
              if node.sequence?
                str += node.sequence
              else
                str += 'short node'
              end
              #log.debug str
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
                  found_direction = true
                end
              else
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
    end
  end
end


require 'bio-logger'
Bio::Log::LoggerPlus.new('finishm')
module Bio
  module FinishM
    module Logging
      def log
        Bio::Log::LoggerPlus['finishm']
      end
    end
  end
end

require 'oligo_designer'
require 'kmer_multi_abundance_file'
require 'kmer_abundance_pattern'

require 'bio-velvet'

require 'assembly/acyclic_connection_finder'
require 'assembly/paths_between_nodes'
require 'assembly/node_finder'
require 'assembly/velvet_graph_sequence_extractor'
require 'assembly/a_b_visualiser'
require 'assembly/coverage_based_graph_filter'
require 'assembly/oriented_node_trail'
require 'assembly/depth_first_search'
require 'assembly/kmer_coverage_based_path_filter'
require 'assembly/graph_generator'
require 'assembly/fluffer'
require 'assembly/scaffold_breaker'
require 'assembly/graph_explorer'
require 'assembly/all_orfs'
require 'assembly/contig_printer'
require 'assembly/bubbly_assembler'
require 'assembly/dijkstra'

require 'finishm/finisher'
require 'finishm/gapfiller'
require 'finishm/wander'
require 'finishm/fluff'
require 'finishm/explore'
require 'finishm/visualise'
require 'finishm/sequence'

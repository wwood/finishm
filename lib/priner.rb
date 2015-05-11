
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
require 'bio-velvet_underground'

require 'assembly/velvet_c_binding'

require 'assembly/acyclic_connection_finder'
require 'assembly/single_coherent_paths_between_nodes'
require 'assembly/node_finder'
require 'assembly/c_probe_node_finder'
require 'assembly/velvet_graph_sequence_extractor'
require 'assembly/hybrid_velvet_graph'
require 'assembly/a_b_visualiser'
require 'assembly/coverage_based_graph_filter'
require 'assembly/oriented_node_trail'
require 'assembly/depth_first_search'
require 'assembly/kmer_coverage_based_path_filter'
require 'assembly/probed_graph'
require 'assembly/read_input'
require 'assembly/graph_generator'
require 'assembly/fluffer'
require 'assembly/scaffold_breaker'
require 'assembly/graph_explorer'
require 'assembly/all_orfs'
require 'assembly/contig_printer'
require 'assembly/dijkstra'
require 'assembly/single_coherent_wanderer'
require 'assembly/connection_interpreter'
require 'assembly/input_genome'
require 'assembly/read_to_node'
require 'assembly/paired_end_neighbour_finder'

require 'assembly/single_ended_assembler'
require 'assembly/bubbly_assembler'
require 'assembly/bad_format_writer'
require 'assembly/height_finder'
require 'assembly/probe_explorer'

require 'finishm/primers'
require 'finishm/primers_check'
require 'finishm/finisher'
require 'finishm/gapfiller'
require 'finishm/wander'
require 'finishm/fluff'
require 'finishm/explore'
require 'finishm/assemble'
require 'finishm/visualise'
require 'finishm/sequence'
require 'finishm/roundup'
require 'finishm/path_counter'

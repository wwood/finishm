
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

require 'kmer_abundance_pattern'

require 'bio-velvet'

require 'assembly/acyclic_connection_finder'
require 'assembly/node_finder'
require 'assembly/velvet_graph_sequence_extractor'
require 'assembly/a_b_visualiser'
require 'assembly/coverage_based_graph_filter'
require 'assembly/oriented_node_trail'


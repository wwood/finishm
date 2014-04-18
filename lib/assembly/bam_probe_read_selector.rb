require 'bio-samtools'

module Bio
  module AssemblyGraphAlgorithms
    class BamProbeReadSelector
      include Bio::FinishM::Logging

      # Given an indexed bam file of reads mapped onto contigs,
      # an array of one or more [contig_name, position, direction] entries (i.e. places in the contigs to locate reads for),
      # a kmer (the match has to be at least one perfect kmer overlapping the position) and a
      # path to a CnyUnifiedSeq.names file, return an Array of read_IDs of reads that can be used to locate the contig
      # ends in the velvet graph.
      #
      # This assumes that velvet hasn't done anything to clean up the graph as cleaning might remove reads
      # of interest
      def find_probes(indexed_bam_file, contig_name_position_directions, kmer, path_to_cny_unified_seq_names_file)
      end

      # Given a contig name and a side, pick out a read that can be used to 'locate'
      # the contig end in the assembly.
      def find_probe_read_name_from_contig_end(indexed_bam_file, contig_name, direction, position, kmer)
        # Search for all reads that overlap the overhang base, and are in the correct direction

        # Reject reads that do not have matching stretches of DNA that are at least kmer length long
        # as these will not be included in the assembly.

        # Return the 'best' read's name and sequence.
      end
    end
  end
end

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
      def find_probes(indexed_bam_file, contig_names_positions_directions, kmer, path_to_cny_unified_seq_names_file)
        # need to check the sequence of the aligned read is the same as what is in the cny_unified_seq_names_file
      end

      # Given a contig name and a side, together with a path to an indexed bam file,
      # pick out a read that can be used to 'locate'
      # the contig end in the assembly, and return a Bio::DB::Alignment object of it
      def find_probe_read_alignment_from_contig_end(indexed_bam_file, contig_name, direction, position, kmer)
        # Search for all reads that overlap the overhang base, and are in the correct direction
        sam = Bio::DB::Sam.new(:bam => indexed_bam_file)
        position_hash = {:chr => contig_name}

        # The probes must overlap the position, to one back from
        # the contig end
        if direction
          position_hash[:start] = position-1
          position_hash[:stop] = position
        else
          position_hash[:start] = position
          position_hash[:stop] = position+1
        end
        sam.each_alignment(position_hash) do |alignment|
          # Reject reads that do not have matching stretches of DNA that are at least kmer length long
          # as these will not be included in the assembly.
          # If it passes, then return the alignment

        end

        # Return the 'best' read's name and sequence.
      end
    end
  end
end

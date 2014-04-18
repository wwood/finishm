module Bio
  module AssemblyGraphAlgorithms
    class BamProbeReadSelector
      include Bio::FinishM::Logging

      # Given a contig name and a side, pick out a read that can be used to 'locate'
      # the contig end in the assembly.
      def find_probe_read_from_contig_end(indexed_bam_file, contig_name, direction, min_overhang, kmer)
        # Search for all reads that overlap the overhang base, and are in the correct direction

        # Reject reads that do not have matching stretches of DNA that are

        # Return the 'best' read.
      end

      #
      def find_probes_from_ends_of_contigs

      end
    end
  end
end

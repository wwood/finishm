require 'bio'

module Bio::AssemblyGraphAlgorithms
  class KmerCoverageBasedPathFilter
    include Bio::FinishM::Logging

    # Remove all paths where the kmer coverage is below the threshold at
    # any point along the path.
    #
    # :paths: an iterable collection of paths
    # :kmer_hash: KmerMultipleAbundanceHash
    # :thresholds: minimum coverages (min numbers of full kmers) required at each point along the path
    #
    # Optional options
    # * :exclude_ending_length => don't filter out paths based on kmers close to the end e.g. :exclude_ending_length => 2 means the first 2 and last 2 kmers are not considered
    def filter(paths, kmer_hash, thresholds, options={})
      # sanity check
      unless kmer_hash.number_of_abundances == thresholds.length
        raise "Unexpectedly found a different number of thresholds and kmer abundance columns"
      end

      passable_paths = []
      paths.each do |path|
        passable = true
        seq = path.sequence

        # remove ends of sequence if kmers don't count at the ends
        if options[:exclude_ending_length]
          ex = options[:exclude_ending_length]
          seq = seq[ex...(seq.length-ex)]
        end

        Bio::Sequence::NA.new(seq).window_search(kmer_hash.kmer_length,1) do |kmer|
          kmer_hash[kmer].each_with_index do |abundance, i|
            if abundance < thresholds[i]
              passable = false
              log.debug "Failing trail #{path.sequence} due to insufficent abundance (#{abundance} from #{kmer_hash[kmer]}) for #{kmer}" if log.debug?
              break
            end
          end
          break if !passable
        end
        passable_paths.push path if passable
      end
      return passable_paths
    end
  end
end

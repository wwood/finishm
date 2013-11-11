require 'csv'

module Bio
  # A class to work with a kmer abundance file format,
  # where the kmer is first, then each abundance after that (no headings, space separated)
  class KmerMultipleAbundanceHash < Hash
    include Bio::FinishM::Logging

    def self.parse_from_file(path)
      obj = self.new
      kmer_length = nil
      num_abundances = nil
      CSV.foreach(path, :col_sep => ' ') do |row|
        kmer = row[0].upcase
        abundances = row[1...row.length]

        kmer_length ||= kmer.length
        if kmer.length != kmer_length
          raise "inconsistent length of kmer found in kmer abundance file, in line: #{row.inspect}"
        end
        num_abundances ||= abundances.length
        if num_abundances != abundances.length
          raise "inconsistent number of abundances found in kmer abundance file, in line: #{row.inspect}"
        end
        obj[kmer] = abundances
      end
    end

    def kmer_length
      each do |kmer, abundances|
        return kmer.length
      end
    end

    def number_of_abundances
      each do |kmer, abundances|
        return abundances.length
      end
    end

    def [](kmer)
      abundances = super(kmer.upcase)
      abundances ||= [0]*number_of_abundances
      return abundances
    end
  end
end

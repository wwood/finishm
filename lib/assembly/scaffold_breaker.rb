require 'bio'

class Bio::FinishM::ScaffoldBreaker
  include Bio::FinishM::Logging

  class UnscaffoldedContig
    attr_accessor :scaffold_position_start, :scaffold_position_end

    # The Scaffold to which this contig once belonged
    attr_accessor :scaffold

    # The actual nucleotide sequence of this contig, from scaffold start position to
    # end (not revcom)
    attr_accessor :sequence

    def length
      @scaffold_position_end - @scaffold_position_start +1
    end

    def name
      contig_number = scaffold.contigs.find_index(self)+1
      if contig_number.nil?
        raise "A contig finds itself unexpectedly not in the scaffold it is supposed to belong to"
      end
      return "#{scaffold.name}_#{contig_number}of#{scaffold.contigs.length}_#{scaffold_position_start}to#{scaffold_position_end}"
    end
  end

  class Scaffold
    # unscaffolded contigs from this scaffold, as an array in sorted order.
    attr_accessor :contigs

    # Name of sequence found in the fasta file
    attr_accessor :name

    # Return an array of Gap objects
    def gaps
      gaps = []
      last_contig = nil
      @contigs.each_with_index do |contig, i|
        if i!=0
          gap = Bio::FinishM::ScaffoldBreaker::Gap.new
          gap.scaffold = self
          gap.start = last_contig.scaffold_position_end + 1
          gap.stop = contig.scaffold_position_start - 1
          gap.number = i-1
          gaps.push gap
        end
        last_contig = contig
      end
      return gaps
    end

    def sequence
      to_return = []
      last_contig = nil
      @contigs.each_with_index do |contig, i|
        if i==0
          to_return.push contig.sequence
        else
          gap_start = last_contig.scaffold_position_end + 1
          gap_stop = contig.scaffold_position_start - 1
          to_return.push 'N'*(gap_stop-gap_start+1)
          to_return.push contig.sequence
        end
        last_contig = contig
      end
      return to_return.join
    end

    # Which contig number is this, in the scaffold?
    def contig_number(contig)
      @contigs.each_with_index do |current_contig, i|
        return i if contig==current_contig
      end
      raise "Contig not found in scaffold"
    end
  end

  class Gap
    attr_accessor :scaffold, :start, :stop, :number

    def coords
      @scaffold.name+':'+(@start+1).to_s+'-'+(@stop).to_s
    end

    #i.e. the number of N characters that would represent this gap
    def length
      @stop-@start+1
    end
  end

  # Given a path to a scaffold fasta file, read in the scaffolds, and break them apart
  # into constituent contigs. Then return an array of Scaffold objects containing the
  # contig information therein.
  def break_scaffolds(contigs_filename)
    scaffolds = []
    Bio::FlatFile.foreach(Bio::FastaFormat, contigs_filename) do |seq|
      scaffold = Scaffold.new
      scaffold.name = seq.definition

      unless seq.seq.match(/^[ATGCN]+$/i)
        example = seq.seq.match(/([^ATGCN])/i)[1]
        log.warn "Found unexpected characters in the sequence #{seq.definition} e.g. #{example}. Replacing them with Ns"
        seq.seq.gsub! /[^ATGCN]/i, 'N'
      end

      if seq.seq.match(/^N+$/i)
        raise "Found a scaffold that contains all N characters, ignoring this (perhaps your input is mangled?): #{scaffold.name}"
      end

      # Find all Ns in the current sequence
      seq.seq.scan(/([^N]+)/i) do
        contig = UnscaffoldedContig.new
        contig.scaffold = scaffold
        contig.scaffold_position_start = $~.offset(0)[0]+1#Convert to 1-based indices in line with bioruby
        contig.scaffold_position_end = $~.offset(0)[1]
        contig.sequence = $~.to_s
        scaffold.contigs ||= []
        scaffold.contigs.push contig
      end
      scaffolds.push scaffold
    end
    log.info "Detected #{scaffolds.length} scaffolds, containing #{scaffolds.collect{|s| s.contigs.length}.reduce(:+)} different contigs"
    return scaffolds
  end
end

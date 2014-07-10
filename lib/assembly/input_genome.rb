class Bio::FinishM::InputGenome
  attr_accessor :scaffolds, :filename, :numbered_probes
  include Bio::FinishM::Logging

  # Return an array of parsed fasta files
  def self.parse_genome_fasta_files(fasta_files, hangover_length, options = {})
    genomes = []
    current_probe_number = 1
    fasta_files.each do |genome_fasta|
      genome = Bio::FinishM::InputGenome.new(
        genome_fasta, hangover_length, :starting_probe_number => current_probe_number
        )
      current_probe_number += genome.number_of_probes

      genomes.push genome
    end
    return genomes
  end

  # Given a fasta file, setup a genome for wandering or gapfilling.
  #
  # Options:
  # :starting_probe_number: number probes starting from this number (default 1)
  def initialize(genome_fasta, hangover_length, options = {})
    starting_probe_number = options[:starting_probe_number]
    starting_probe_number ||= 1

    @filename = genome_fasta
    scaffolds = Bio::FinishM::ScaffoldBreaker.new.break_scaffolds(genome_fasta)
    @scaffolds = remove_overly_short_contigs!(scaffolds, hangover_length)

    # Remove scaffolds that have no good
    num_too_short_scaffolds = 0
    @scaffolds.reject! do |scaff|
      rej = scaff.contigs.empty?
      num_too_short_scaffolds += 1 if rej
      rej
    end
    if num_too_short_scaffolds > 0
      log.warn "Removed #{num_too_short_scaffolds} scaffolds entirely as they were too short (or made up of all short contigs)"
    end

    generate_numbered_probes(hangover_length, starting_probe_number)
  end

  def remove_overly_short_contigs!(scaffolds, hangover_length)
    num_contigs_removed = 0
    scaffolds.each do |scaffold|
      num_original_contigs = scaffold.contigs.length
      scaffold.contigs.reject! do |contig|
        contig.length < 2*hangover_length
      end
      num_contigs_removed += num_original_contigs - scaffold.contigs.length
    end
    if num_contigs_removed > 0
      log.warn "Removed #{num_contigs_removed} contigs from within scaffolds as they were too short"
    end
    return scaffolds
  end

  def generate_numbered_probes(overhang, starting_probe_number)
    @numbered_probes = []
    @probe_number_to_scaffold_and_contig_and_side = {}

    current_probe_number = starting_probe_number
    overly_short_sequence_count = 0
    @scaffolds.each_with_index do |scaffold, scaffold_index|
      scaffold.contigs.each_with_index do |contig, contig_index|
        if contig.sequence.length < 2*overhang
          log.warn "Not attempting to make connections from overly short contig: it is the #{contig_index+1}th contig in scaffold `#{scaffold.name}' from the genome in `#{@filename}')"
          overly_short_sequence_count += 1
          nil
        else
          sequence = contig.sequence

          probe1 = NumberedProbe.new
          probe1.contig = contig
          probe1.number = current_probe_number; current_probe_number += 1
          probe1.side = :start
          fwd2 = Bio::Sequence::NA.new(sequence[0...overhang])
          probe1.sequence = fwd2.reverse_complement.to_s

          probe2 = NumberedProbe.new
          probe2.contig = contig
          probe2.number = current_probe_number; current_probe_number += 1
          probe2.side = :end
          probe2.sequence = sequence[(sequence.length-overhang)...sequence.length].to_s

          @numbered_probes[scaffold_index] ||= []
          @numbered_probes[scaffold_index][contig_index] = [probe1, probe2]

          @probe_number_to_scaffold_and_contig_and_side[probe1.number] = [scaffold, contig, :start]
          @probe_number_to_scaffold_and_contig_and_side[probe2.number] = [scaffold, contig, :end]
        end
      end
    end
    log.debug "Generated #{current_probe_number-starting_probe_number} probes for #{@filename}" if log.debug?
    if overly_short_sequence_count > 0
      log.warn "Skipping #{overly_short_sequence_count} overly short contigs" if log.warn?
    end
  end

  def number_of_probes
    @numbered_probes.flatten.length
  end

  def each_numbered_probe
    @numbered_probes.flatten.each do |probe|
      yield probe
    end
  end

  def each_gap_probe_pair(scaffold_index)
    last_probe_pair = nil
    @numbered_probes[scaffold_index].each do |probe_pair|
      unless probe_pair.nil?
        unless last_probe_pair.nil?
          yield last_probe_pair[1], probe_pair[0]
        end
        last_probe_pair = probe_pair
      end
    end
  end

  def each_scaffold_end_numbered_probe
    @numbered_probes.each_with_index do |scaffold_indices, i|
      # yield the first and last probe of this scaffold
      yield scaffold_indices[0][0]
      yield scaffold_indices[-1][1]
    end
  end

  def first_probe(scaffold_index)
    @numbered_probes[scaffold_index][0][0]
  end

  def last_probe(scaffold_index)
    @numbered_probes[scaffold_index][-1][1]
  end

  # Return true if probe number given is the probe at the beginning of the scaffold
  # or false if it is at the end. raise if unknown.
  def probe_at_start_of_scaffold?(probe_index)
    scaffold, contig, side = @probe_number_to_scaffold_and_contig_and_side[probe_index]
    if side == :start
      return true
    elsif side == :end
      return false
    else
      raise
    end
  end

  def probe_sequences
    seqs = []
    each_numbered_probe do |probe|
      unless probe.nil?
        seqs.push probe.sequence
      end
    end
    return seqs
  end

  # The length of the gap between contig i and contig i+1 from the specified scaffold
  # (both numbers are 0-based indices)
  def gap_length(scaffold_index, gap_index)
    @scaffolds[scaffold_index].gaps[gap_index].length
  end
end

class NumberedProbe
  attr_accessor :number, :contig, :side, :sequence

  def index
    @number - 1
  end
end

class OligoDesigner
  # Given a sequence, find the subsequence that starts at the 5' end (ie the
  # start of the string), and ends when the melting temperature is maximal but
  # below the max_temperature requires the oligotm program to be available on the
  # cmd line.
  #
  # * nucleotide_string: the full sequence that we are choosing oligos from
  # * max_temperature: the maximal temperature to start things off at
  # * gc_clamp: require this many G or C residues at the 3' end of the oligo.
  def just_below(nucleotide_string, max_temperature, gc_clamp=0)
    # initial conditions
    guess = 0
    guess_temp = 0

    # loop around
    while guess_temp < max_temperature
      guess += 1
      guess_temp = melting_temperature nucleotide_string[0..guess-1]

      # break if there's we've reached the end of the line
      return nucleotide_string if guess > nucleotide_string.length
    end
    return nucleotide_string[0..guess-2]
  end

  # Rank oligomers within some constraints.
  def possible_oligos_ordered_by_temperature_difference(nucleotide_string, min_temperature, best_temperature, max_temperature, gc_clamp)
    default_distance = lambda do |seq, tm|
    # fails constraints if not enough GC clamp
      if seq[seq.length-gc_clamp..seq.length-1].gsub(/[gc]/i,'').length > 0
      false
      else
      # the sequence is within contraints. The melting temperature closest to the best wins.
      tm_diff = (best_temperature-tm).abs
      tm_diff
      end
    end

    # initial conditions
    guess = 0
    guess_temp = 0
    # arrays to fill with possible possibles
    oligos = []

    # loop around, until max temperature is reached
    while guess_temp < max_temperature
      guess += 1
      seq = nucleotide_string[0..guess-1]
      guess_temp = melting_temperature seq

      # Add it to the list if there is enough temperature
      if guess_temp > min_temperature and guess_temp < max_temperature
      o = Oligo.new
      o.sequence = seq
      o.tm = guess_temp
      oligos.push o
      end

      # break if there's we've reached the end of the line
      break if guess > nucleotide_string.length-1
    end

    # Convert sequences into distances
    oligos.each do |oligo|
      oligo.distance = default_distance.call(oligo.sequence, oligo.tm)
    end

    # remove sequences that don't meet the constraints, and sort the rest with
    # smallest distance first
    return oligos.reject{|o| o.distance == false}.sort{|a,b|
      a.distance<=>b.distance
    }.collect{|o| o.sequence}
  end
  alias_method :order, :possible_oligos_ordered_by_temperature_difference

  # A simple method to return the melting temperature of a particular nucleotide
  # string. Uses oligotm on the command line.
  def melting_temperature(nucleotide_string)
    #`oligotm -tp 1 -sc 1 -n 0.8 -d 500 -mv 0 -dv 50 '#{nucleotide_string}'`.to_f
    `oligotm -tp 1 -sc 1 -n 0.2 -d 2 -mv 1 '#{nucleotide_string}'`.to_f
  end

  private

  class Oligo
    attr_accessor :sequence, :tm, :distance
  end
end
# A pattern of presence/absence/neither across a run of kmers
class KmerAbundancePattern < Array
  def binary_string
    to_return = ''
    each do |present|
      to_return += case present
      when true
        '1'
      when false
        '0'
      when '-'
        '-'
      else
        raise "Unexpected pattern atom found: #{present}"
      end
    end
    to_return
  end

  # Parse a 100001011 type representation
  def parse_from_human(boolean_pattern)
    boolean_pattern.each_char do |char|
      if char == '1'
        push true
      elsif char == '0'
        push false
      else
        raise "Unexpected pattern character: #{char}"
      end
    end
  end

  def consistent_with(another_pattern)
    unless length == another_pattern.length
      raise "Unexpected comparison of this pattern #{inspect} with another: #{another_pattern.inspect}"
    end
    each_with_index do |bool, i|
      return false if bool != another_pattern[i]
    end
    return true
  end

  def parse_from_kmer_abundance(abundances, lower_limit, upper_limit)
    abundances.each do |a|
      if a>=upper_limit
        push true
      elsif a<=lower_limit
        push false
      else
        push '-'
      end
    end
  end
end

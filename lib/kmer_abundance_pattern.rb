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
    self[0...length] = [] #remove the last pattern if it existed
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

  # Return true if this pattern is exactly the same
  # as another pattern
  #
  # e.g. 101 is same_as? 101 but not 111 or 110
  def same_as?(another_pattern)
    unless length == another_pattern.length
      raise "Unexpected comparison of this pattern #{inspect} with another: #{another_pattern.inspect}"
    end
    each_with_index do |bool, i|
      return false if bool != another_pattern[i]
    end
    return true
  end

  # Return true if another_pattern shows presence in all places
  # where this pattern is present, (but maybe more)
  #
  # e.g. 101 is consisten with 101 and 111, but not 011
  #
  # Behaviour not defined when in no-mans land
  def consistent_with?(another_pattern)
    unless length == another_pattern.length
      raise "Unexpected comparison of this pattern #{inspect} with another: #{another_pattern.inspect}"
    end
    each_with_index do |bool, i|
      return false if bool and !another_pattern[i]
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

class Bio::AssemblyGraphAlgorithms::PairedEndAssembler < Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  # Assemble considering reads as a possibly paired-ended. Options are as per SingleEndedAssembler#assemble_from,
  # with the addition of
  # :min_insert_size: minimum length of fragment pair required to satisfy the additional
  # constraints of the paired end assembler.
  # :max_insert_size: maximum length of fragment pair.
  def assemble_from(initial_path, visited_nodes)
    visited_nodes = Set.new
    while true
      # Try to assemble using single ended techniques first, and only if that fails fall
      # back to paired-end techniques.
      path, visited_nodes, next_neighbours = super(initial_path, visited_nodes)

      # The next_neighbours
      if next_neighbours.empty?
        # No-where to go, do nothing
        # TODO: try to jump over the gap using paired-end sequences
      elsif next_neighbours.length < 2
        raise "Programming error"
      else
        # Choose between forks based on paired-end data
        next_neighbours.select! do |oneigh|
          confirm_connection_backwards(oneigh, path)
        end
      end
    end
  end

  # Return true if the backward connection is strong enough to warrant adding
  # the new_node to the current_path. Else return false
  #
  # In order to qualify as a warranted path, reads from the pair must have
  # at least 1 connection backwards at least options[:min_insert_size] backwards
  # from the end of the current path, but not more than options[:max_insert_size]
  # Connection length is the length of the read pair's insert size.
  def confirm_connection_backwards(new_onode, current_path, min_insert_size, max_insert_size)
    min_insert_size = @assembly_options[:min_insert_size]
    max_insert_size = @assembly_options[:max_insert_size]

    # Collect the reads that would qualify if they were in the new node
    qualifying_reads = Set.new
    current_path.reverse.each do |node|

    end

    # Look at all the sequences in the new onode. Do any of them qualify
    new_onode.node.short_reads.each do |short|

    end
  end
end

class Bio::FinishM::ORFsFinder
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :min_orf_length => 100
    }

  def add_options(optparse_object, options)
    options.merge! Bio::FinishM::Visualise::DEFAULT_OPTIONS
    options.merge! DEFAULT_OPTIONS
    optparse_object.banner = "\nUsage: finishm find_orfs --assembly-???

    Find possible open reading frames in assembly graph
    \n\n"

    optparse_object.separator "Input genome information"
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-exploration arguments:\n\n"
    Bio::FinishM::Visualise.new.add_probe_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options(optparse_object, options)
  end

  def validate_options(options, argv)
    visualise = Bio::FinishM::Visualise.new
    return visualise.validate_argv_length(argv) ||
      visualise.validate_probe_options(options) ||
      visualise.validate_assembly_options(options)
  end

  def run(options, argv)
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    visualise = Bio::FinishM::Visualise.new

    if options[:interesting_probes] or options[:interesting_probe_names]
      finishm_graph, interesting_node_ids = visualise.generate_graph_from_probes(read_input, options)
    elsif options[:interesting_nodes]
      finishm_graph = visualise.generate_graph_from_nodes(read_input, options)
      interesting_node_ids = options[:interesting_nodes]
    elsif options[:assembly_files]
      finishm_graph, interesting_node_ids, = visualise.generate_graph_from_assembly(read_input, options)
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end

    if options[:graph_search_leash_length]
      #log.info "Finding nodes within the leash length of #{options[:graph_search_leash_length] }.."
      nodes_within_leash, node_ids_at_leash = visualise.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
      log.info "Found #{node_ids_at_leash.length} nodes at the end of the #{options[:graph_search_leash_length] }bp leash" if options[:graph_search_leash_length]

      options[:range] = nodes_within_leash
    else
      options[:range] = finishm_graph.graph.nodes
    end

    initial_onodes = Bio::FinishM::PathCounter.new.get_leash_start_nodes(finishm_graph, options[:range])
    find_orfs_in_graph(finishm_graph, initial_onodes, options)
  end

  def find_orfs_in_graph(finishm_graph, initial_onodes, options={})
    initial_paths = initial_onodes.collect do |onode|
      path = Bio::Velvet::Graph::OrientedNodeTrail.new
      path.add_oriented_node onode
      path
    end

    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    orf_trails = orfer.find_orfs_in_graph(finishm_graph.graph, initial_paths,
        options[:min_orf_length], options[:range])

    found_orfs = {}
    orf_trails.each do |path|
      sequences, sequences = path.otrail.sequence_within_path
      (0..1).each do |dir|
        if dir == 0
          result = path.fwd_orfs_result
          trail = path.otrail.trail
          sequence = fwd_trail
        else
          result = path.twin_orfs_result
          trail = path.otrail.reverse.trail
        end

        result.start_stop_pairs.each do |pair|
          start_index, length_to_start, end_index, length_to_end = orf_end_nodes(trail, pair[0], pair[1])
          start_position_in_trail = pair[0] - length_to_start
          end_position_in_trail = pair[1] - length_to_start
          orf_trail = trail[start_index..end_index]
          orf_name = "#{orf_trail.collect{|onode| ondoe.to_shorthand}.join(',')}[#{start_position_in_trail}:#{end_position_in_trail}]"
          next if found_orfs.has_key? orf_name
          orf_sequence = sequences[dir][pair[0]-2, pair[1]]
          found_orfs[orf_name] = orf_sequence
        end
      end
    end

    found_orfs.each_pair do |name, sequence|
      puts ">#{name}"
      puts sequence
      puts
    end
  end

  def orf_end_nodes(path, start_position, stop_position)
    start_of_current_index = 0
    index = 0

    indices = [first, last].collect do |position|
      while current_node = path[index]
        start_of_next_index = start_of_current_index + current_node.node.length_alone
        if position_of_next < position
          index += 1
          start_of_current_index = start_of_next_index
          next
        end

        return [index, start_of_current_index]
      end
      raise "Position #{position} not in path" if current_node.nil?
    end

    return indices.flatten
  end

  def orf_to_settable(path, start_index, start_offset, end_index, end_offset)
    [path[start_index..end_index].collect{|onode| onode.to_settable},[start_offset, end_offset]]
  end
end

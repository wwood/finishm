class Bio::FinishM::ORFsFinder
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :min_orf_length => 100
    }

  def add_options(optparse_object, options)
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

    if options[:interesting_probes] or options[:interesting_probe_names]
      finishm_graph, interesting_node_ids = generate_graph_from_probes(read_input, options)
    elsif options[:interesting_nodes]
      finishm_graph = generate_graph_from_nodes(read_input, options)
      interesting_node_ids = options[:interesting_nodes]
    elsif options[:assembly_files]
      finishm_graph, interesting_node_ids, = generate_graph_from_assembly(read_input, options)
    else
      finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)
    end

    initial_path = Bio::Velvet::Graph::OrientedNodeTrail.create_from_shorthand options[:start_node], finishm_graph

    if options[:end_nodes]
      options[:terminal_nodes] = options[:end_nodes].collect{|id| finishm_graph.nodes[id]}
    end

    log.info "Finding nodes from which to begin search.."
    nodes_within_leash, node_ids_at_leash = visualise.get_nodes_within_leash(finishm_graph, interesting_node_ids, options)
    log.info "Found #{start_node_ids.length} nodes"



    find_orfs_in_graph(finishm_graph, initial_path, options)
  end

  def get_start_and_end_nodes_within_leash(finishm_graph, node_ids, options)
    log.info "Finding nodes within the leash length of #{options[:graph_search_leash_length] } with maximum node count #{options[:max_nodes] }.."
    dijkstra = Bio::AssemblyGraphAlgorithms::Dijkstra.new

    @finder = Bio::FinishM::PairedEndNeighbourFinder.new(finishm_graph, 500) #TODO: this hard-coded 100 isn't great here
    @finder.min_adjoining_reads = options[:min_adjoining_reads]
    @finder.max_adjoining_node_coverage = options[:max_adjoining_node_coverage]

    nodes_within_leash_hash = dijkstra.min_distances_from_many_nodes_in_both_directions(
      finishm_graph.graph, node_ids.collect{|n| finishm_graph.graph.nodes[n]}, {
        :ignore_directions => true,
        :leash_length => options[:graph_search_leash_length],
        :max_nodes => options[:max_nodes],
        :neighbour_finder => @finder
        })
    nodes_within_leash = nodes_within_leash_hash.keys.collect{|k| finishm_graph.graph.nodes[k[0]]}
    log.info "Found #{nodes_within_leash.collect{|o| o.node_id}.uniq.length} node(s) within the leash length"

    # These nodes are at the end of the leash - a node is in here iff
    # it has a neighbour that is not in the nodes_within_leash
    start_node_ids = Set.new
    end_node_ids = Set.new
    nodes_within_leash_hash.keys.each do |node_and_direction|
      # Add it to the set if 1 or more nieghbours are not in the original set
      node = finishm_graph.graph.nodes[node_and_direction[0]]
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new node, node_and_direction[1]
      onode.next_neighbours(finishm_graph.graph).each do |oneigh|
        if !nodes_within_leash_hash.key?(oneigh.to_settable)
          node_ids_at_leash << node_and_direction
          break #it only takes one to be listed
        end
      end
    end


    return nodes_within_leash.uniq, start_node_ids, end_node_ids
  end

  def find_orfs_in_graph(finishm_graph, initial_path, options={})
    orfer = Bio::AssemblyGraphAlgorithms::AllOrfsFinder.new
    orf_trails = orfer.find_orfs_in_connection(finishm_graph.graph, initial_path, options[:terminal_nodes],
        options[:min_orf_length], options[:leash_length], options[:max_gapfill_paths], options[:max_cycles],
        options[:max_explore_nodes])

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

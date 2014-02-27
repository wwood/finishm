class Bio::FinishM::Sequence
  include Bio::FinishM::Logging

  class PathSteppingStone
    attr_accessor :node_id, :first_side
  end

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm sequence --assembly-??? --path PATH

    Given a series of nodes and orientations, print the DNA sequence of the given path
    \n\n"

    options.merge!({
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--path PATH", Array, "A comma separated list of node IDs and orientations - explore from these probe IDs in the graph e.g. '4s,2s,3e' means start at the start of node 4, connecting to the beginning of node 2 and finally the end of probe 3. [required]") do |arg|
      options[:path] = arg.collect do |str|
        if matches = str.match(/^([01-9]+)([se])$/)
          stone = PathSteppingStone.new
          stone.node_id = matches[1].to_i
          if matches[2] == 's'
            stone.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
          elsif matches[2] == 'e'
            stone.first_side = Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
          else
            raise "programming error"
          end
          stone
        else
          raise "Unable to parse stepping stone along the path: `#{arg}'. Entire path was `#{arg}'."
        end
      end
    end
    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional graph-related arguments:\n\n"
    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      if options[:path].nil?
        return "No path defined, so don't know how to procede through the graph"
      end

      # Need reads unless there is already an assembly
      unless options[:previous_assembly] or options[:previously_serialized_parsed_graph_file]
        return Bio::FinishM::ReadInput.new.validate_options(options, [])
      else
        return nil
      end
    end
  end

  def run(options, argv)
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options

    # Generate the assembly graph
    log.info "Reading in or generating the assembly graph"
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph([], read_input, options)

    # Build the oriented node trail
    log.info "Building the trail from the nodes"
    trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    options[:path].each do |stone|
      log.debug "Adding stone to the trail: #{stone.inspect}"
      node = finishm_graph.graph.nodes[stone.node_id]
      if node.nil?
        raise "Unable to find node ID #{stone.node_id} in the graph, so cannot continue"
      end

      # check that the path actually connects in the graph, otherwise stop.
      is_neighbour = false
      unless trail.length == 0 #don't worry about the first stepping stone
        trail.neighbours_of_last_node(finishm_graph.graph).each do |oneigh|
          log.debug "Considering neighbour #{oneigh.inspect}"
          is_neighbour = true if oneigh.node == node and oneigh.first_side == stone.first_side
        end
        unless is_neighbour
          raise "In the graph, the node #{trail.last.to_s} does not connect with #{stone.inspect}"
        end
      end

      # OK, all the checking done. Actually add it to the trail
      trail.add_node node, stone.first_side
    end

    # Print the sequence
    puts trail.sequence
  end
end

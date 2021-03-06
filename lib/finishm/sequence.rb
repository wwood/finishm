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

    # Parse a string like '4s,2s,3e' into a programmatic version of a path
    parse_path_string = lambda do |path_string|
      path_string.collect do |str|
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

    optparse_object.separator "\nOne of the following path defining arguments must be defined:\n\n"
    optparse_object.on("--path-ids PATH", "A comma separated list of node IDs - the program attempts to determine the orientations automatically") do |arg|
      options[:path_ids] = arg
    end
    optparse_object.on("--path PATH", Array, "A comma separated list of node IDs and orientations - explore from these probe IDs in the graph e.g. '4s,2s,3e' means start at the start of node 4, connecting to the beginning of node 2 and finally the end of probe 3.") do |arg|
      options[:paths] = [parse_path_string.call(arg)]
    end
    optparse_object.on("--paths PATHS", "A colon separated list of comma separated lists of node IDs and orientations - e.g. '4s,2s,3e:532s,465s' means print 2 different paths") do |arg|
      raise "Only one of --path and --paths can be specified" unless options[:paths].nil?
      options[:paths] = []
      arg.split(':').each do |split|
        split.strip!
        next if split == ''
        options[:paths].push parse_path_string.call(split.split(','))
      end
      log.info "Read in #{options[:paths] } path definitions"
      if log.debug?

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
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      if options[:path_ids]
        if options[:paths]
          return "Multiple ways to define the path given, one at a time please"
        end
      else
        if options[:paths].nil? or options[:paths].empty?
          return "No path defined, so don't know how to procede through the graph"
        end
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

    print_trail = lambda do |oriented_trail|
      print '>'
      puts oriented_trail.to_shorthand
      puts oriented_trail.sequence
    end

    if options[:path_ids]
      trail = Bio::Velvet::Graph::OrientedNodeTrail.create_from_super_shorthand(options[:path_ids], finishm_graph.graph)
      print_trail.call trail

    else
      # Build the oriented node trail
      log.info "Building the trail(s) from the nodes"
      options[:paths].each do |path|
        trail = Bio::Velvet::Graph::OrientedNodeTrail.new
        path.each do |stone|
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
        print_trail.call trail
      end
    end
  end
end

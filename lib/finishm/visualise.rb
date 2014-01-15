class Bio::FinishM::Visualiser
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm visualise --assembly-??? <output_visualisation_file>

    Visualise an assembly graph
    \n\n"

    options.merge!({
    })
    optparse_object.separator "Output visualisation formats (one or more of these must be used)"
    optparse_object.on("--assembly-png PATH", "Output assembly as a PNG file [default: off]") do |arg|
      options[:output_graph_png] = arg
    end
    optparse_object.on("--assembly-svg PATH", "Output assembly as a SVG file [default: off]") do |arg|
      options[:output_graph_svg] = arg
    end
    optparse_object.on("--assembly-dot PATH", "Output assembly as a DOT file [default: off]") do |arg|
      options[:output_graph_dot] = arg
    end

    optparse_object.separator "\nIf an assembly is to be done, there must be some definition of reads:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
#    optparse_object.on("--contigs FILE", "Fasta file containing contigs to find the fluff on [required]") do |arg|
#      options[:contigs_file] = arg
#    end
#    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
#      options[:contig_end_length] = arg.to_i
#    end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      if options[:output_graph_png].nil? and options[:output_graph_svg].nil? and options[:output_graph_dot]
        return "No visualisation output format/file given, don't know how to visualise"
      end
      #TODO: this needs to be improved.

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

    # Display the output graph visualisation
    log.info "Converting assembly to a graphviz"
    viser = Bio::Assembly::ABVisualiser.new
    gv = viser.graphviz(finishm_graph.graph, {:start_node_ids => finishm_graph.probe_nodes.collect{|node| node.node_id}})

    if options[:output_graph_png]
      log.info "Writing PNG #{options[:output_graph_png]}"
      gv.output :png => options[:output_graph_png], :use => :neato
    end
    if options[:output_graph_svg]
      log.info "Writing SVG #{options[:output_graph_svg]}"
      gv.output :svg => options[:output_graph_svg], :use => :neato
    end
    if options[:output_graph_dot]
      log.info "Writing DOT #{options[:output_graph_dot]}"
      gv.output :dot => options[:output_graph_dot] if options[:output_graph_dot]
    end
  end
end

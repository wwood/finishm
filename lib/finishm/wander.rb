class Bio::FinishM::Wanderer
  include Bio::FinishM::Logging

  # Collect desciptions about the probes so that they can be inspected more easily given a probe index
  class ProbeDescription
    attr_accessor :sequence_name, :side

    def to_s
      "#{@sequence_name}.#{@side}"
    end
  end

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm wander --contigs <contig_file> --fastq-gz <reads..> --output-connections <output.csv>

    Takes a set of contigs/scaffolds and finds connections in the graph between them. A connection here is given as
    the length of the shortest path between them, without actually computing all the paths.

    This method can be used for 'pre-scaffolding', in the following sense. If the shortest path between
    two contig ends is 10kb, and a mate pair library with insert size 2kb suggests a linkage
    between the two ends, then the mate pair linkage is likely false (as long as there is sufficient
    coverage in the reads, and not overwhelming amounts of strain heterogeneity, etc.).

    \n\n"

    options.merge!({
      :contig_end_length => 200,
      :graph_search_leash_length => 20000,
      :unscaffold_first => false,
    })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
      options[:contigs_file] = arg
    end
    optparse_object.on("--output-connections PATH", "Output connections in tab-separated format [required]") do |arg|
      options[:output_connection_file] = arg
    end
    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length]}]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: #{options[:graph_search_leash_length]}]") do |arg|
      options[:graph_search_leash_length] = arg
    end
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs[default: #{options[:unscaffold_first]}]") do |arg|
      options[:unscaffold_first] = true
    end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0]}"
    else
      [
        :contigs_file,
        :output_connection_file
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv)
    # Read in all the contigs sequences, removing those that are too short
    probe_sequences = []
    sequence_names = []
    process_sequence = lambda do |name, seq|
      if seq.length < 2*options[:contig_end_length]
        log.warn "Not attempting to make connections from this contig, as it is overly short: #{name}"
        nil
      else
        sequence_names.push name

        sequence = seq.seq
        fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
        probe_sequences.push fwd2.reverse_complement.to_s

        probe_sequences.push sequence[(sequence.length-options[:contig_end_length])...sequence.length]

        # 'return' the probe indices that have been assigned
        [probe_sequences.length-2, probe_sequences.length-1]
      end
    end

    scaffolds = nil #Array of Bio::FinishM::ScaffoldBreaker::Scaffold objects. Only set when options[:unscaffold_first] is true
    scaffolded_contig_to_probe_ids = {}
    if options[:unscaffold_first]
      log.info "Unscaffolding scaffolds (before trying to connect them together again)"
      scaffolds = Bio::FinishM::ScaffoldBreaker.new.break_scaffolds options[:contigs_file]
      scaffolds.each do |scaffold|
        scaffold.contigs.each do |contig|
          process_sequence.call contig.name, contig.sequence
        end
      end
    else
      # Else don't split up any of the sequences
      log.info "Reading input sequences.."
      Bio::FlatFile.foreach(options[:contigs_file]) do |seq|
        process_sequence.call seq.definition, seq.seq
      end
    end
    log.info "Searching from #{probe_sequences.length} different contig ends (#{probe_sequences.length / 2} contigs)"

    # Generate the graph with the probe sequences in it.
    read_input = Bio::FinishM::ReadInput.new
    read_input.parse_options options
    finishm_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    # Loop over the ends, trying to make connections from each one
    cartographer = Bio::AssemblyGraphAlgorithms::AcyclicConnectionFinder.new

    log.info "Finding possible connections with a depth first search"
    first_connections = cartographer.depth_first_search_with_leash(finishm_graph, options[:graph_search_leash_length])
    log.info "Found #{first_connections.length} connections with less distance than the leash length, out of a possible #{probe_sequences.length*(probe_sequences.length-1) / 2}"

    probe_descriptions = []
    (0...finishm_graph.probe_nodes.length).each do |i|
      desc = ProbeDescription.new
      if i % 2 == 0
        desc.side = 'start'
        desc.sequence_name = sequence_names[i / 2]
      else
        desc.side = 'end'
        desc.sequence_name = sequence_names[(i-1) / 2]
      end
      probe_descriptions.push desc
    end

    # Write out connections to the given file
    File.open(options[:output_connection_file], 'w') do |out|
      first_connections.each do |node_indices, distance|
        sequence_names_and_directions = node_indices.collect do |i|
          probe_descriptions[i].to_s
        end
        out.puts [
          sequence_names_and_directions,
          distance
        ].flatten.join("\t")
      end
    end

    # Make an undirected graph to represent the connections so it is easier to work with
#     graph = Yargraph::UndirectedGraph.new
#     all_probe_indices = (0...finishm_graph.probe_nodes.length).to_a
#     all_probe_indices.each{|i| graph.add_vertex i}
#     first_connections.keys.each{|join| graph.add_edge join[0], join[1]}

#     # Print out which contig ends have no connections
#     unconnected_probe_indices = []
#     graph.vertices.each do |vertex_i|
#       unconnected_probe_indices.push vertex_i if graph.degree(vertex_i) == 0
#     end
#     log.info "Found #{unconnected_probe_indices.length} contig ends that had no connection to any other contig"
#     unless unconnected_probe_indices.empty?
#       log.warn "Unconnected contig ends such as this indicate an error with the assembly - perhaps it is incomplete or there has been a misassembly."
#     end
#     unconnected_probe_indices.each do |unconnected|
#       log.warn "#{probe_descriptions[i]} was not connected to any other contig ends"
#     end

#     # find probes with exactly one connection
#     singly_connected_indices = []
#     graph.vertices.each do |vertex_i|
#       singly_connected_indices.push vertex_i if graph.degree(vertex_i) == 1
#     end
#     log.info "Found #{singly_connected_indices.length} contig ends that connect to exactly one other contig ends, likely indicating that they can be scaffolded together:"
#     singly_connected_indices.each do |i|
#       neighbour = graph.neighbours[i][0] #there must be only 1 neighbour by the definition of degree
#       desc = probe_descriptions[i]
#       neighbour_desc = probe_descriptions[neighbour]
#       log.info "The connection between #{desc} and #{neighbour} appears to be a good one"
#     end

#     # If we are working with a scaffold, compare the original scaffolding with graph
#     # as it now is
#     if options[:unscaffold_first]
#       # Of each connection in the scaffold, is that also an edge here? One would expect so given a sensible leash length
#       scaffolds.each do |scaffold|
#         last_contig = nil
#         scaffolds.contigs.each_with_index do |contig, contig_index|
#         end
#       end
#     end

#     #TODO: implemented this in repo hamiltonian_cycler, need to incorporate it here. See also a script in luca/bbbin that uses that library.
#     #TODO: look for hamiltonian paths as well as hamiltonian cycles

  end
end

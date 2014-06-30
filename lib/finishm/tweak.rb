class Bio::FinishM::Wanderer
  include Bio::FinishM::Logging

  DEFAULT_OPTIONS = {
    :contig_end_length => 200,
    :graph_search_leash_length => 20000,
    :unscaffold_first => false,
    :recoherence_kmer => 1,
    }

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm roundup --genomes <genome1.fasta>[,<genome2.fasta>,...] --fastq-gz <reads..> --output-directory <directory>

    Takes one or more genomes and tries to improve their quality by reducing the number of
    scaffolds and N characters they contain.

    Example:

      finishm roundup --genomes genome1.fasta,genome2.fasta --fastq-gz reads.1.fq.gz,reads.2.fq.gz --output-directory finishm_roundup_results

    That will create a collapsed de-Bruijn graph from reads.1.fq.gz and reads.2.fq.gz, then try to find connections between
    the starts and the ends of the contigs in genome1.fasta. If any connections between contigs are mutually exclusive,
    then they are incorporated into scaffolds together, and gapfilling is attempted. The final sequences are output in
    the finishm_roundup_results directory in FASTA format. The procedure is then repeated for genome2.fasta.

    \n\n"

    options.merge!(DEFAULT_OPTIONS)

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--contigs FILE", Array, "fasta file of single contig containing Ns that are to be closed [required]") do |arg|
      options[:assembly_files] = arg
    end
    optparse_object.on("--output-directory PATH", "Output results to this directory [required]") do |arg|
      options[:output_directory] = arg
    end

    optparse_object.separator "\nThere must be some definition of reads too:\n\n" #TODO improve this help
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--overhang NUM", Integer, "Start assembling this far from the ends of the contigs [default: #{options[:contig_end_length] }]") do |arg|
      options[:contig_end_length] = arg.to_i
    end
    optparse_object.on("--recoherence-kmer NUM", Integer, "Use a kmer longer than the original velvet one, to help remove bubbles and circular paths [default: none]") do |arg|
      options[:recoherence_kmer] = arg
    end
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this far and not much more [default: #{options[:graph_search_leash_length] }]") do |arg|
      options[:graph_search_leash_length] = arg
    end
    optparse_object.on("--unscaffold-first", "Break the scaffolds in the contigs file apart, and then wander between the resultant contigs. This option is only relevant to the wander step; gapfilling is attempted on all gaps by default. [default: #{options[:unscaffold_first] }]") do
      options[:unscaffold_first] = true
    end
    #optparse_object.on("--proceed-on-short-contigs", "By default, when overly short contigs are encountered, finishm croaks. This option stops the croaking [default: #{options[:proceed_on_short_contigs] }]") do
    #  options[:proceed_on_short_contigs] = true
    #end

    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      [
        :assembly_files,
        :output_directory,
      ].each do |sym|
        if options[sym].nil?
          return "No option found to specify #{sym}."
        end
      end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv=[])
    # Make sure output directory is writeable to avoid late croaking
    output_directory = setup_output_directory options[:output_directory]

    # Gather the probes from each genome supplied
    genomes = []
    current_probe_number
    options[:assembly_files].each do |genome_fasta|
      genome = Genome.new
      genome.filename = genome_fasta
      genome.scaffolds = Bio::FinishM::ScaffoldBreaker.break_scaffolds(genome_fasta)

      genome.generate_numbered_probes(options[:contig_end_length], current_probe_number)
      current_probe_number += genome.number_of_probes

      genomes.push genome
    end

    # Generate one velvet assembly to rule them all (forging the assembly is hard work..)
    probe_sequences = genomes.collect{|genome| genome.numbered_probes}.flatten.collect{|probe| probe.sequence}
    # Generate the graph with the probe sequences in it.
    read_input.parse_options options
    master_graph = Bio::FinishM::GraphGenerator.new.generate_graph(probe_sequences, read_input, options)

    gapfill_scaffold = lambda do |scaffold|
      connections = []
      probe_pairs = scaffold.gap_probe_pairs
      probe_pairs.each do |probe1, probe2|
        connections.push gapfiller.gapfill(master_graph, probe1.number, probe2.number)
      end
      connections
    end

    # For each genome, wander, gapfill, then output
    wanderer = Bio::FinishM::Wanderer.new
    gapfiller = Bio::FinishM::GapFiller.new
    genomes.each do |genome|
      # wander using just the probes on the ends of the scaffolds
      connected_scaffolds = wander_a_genome(wanderer, genome, master_graph, options)

      # gapfill between
      # (1) interpreted_connections
      # (2) gaps that were present before above wander
      connections = []
      connected_scaffolds.each do |cross_scaffold_connection|
        # Gapfill contigs within the scaffold on the extreme LHS
        connections.push gapfill_scaffold.call(scaffold)

        scaffold.gaps.each_with_index do |gap, i|
          # Gapfill across the new gap between scaffolds
          gapfiller.gapfill(master_graph, scaff)

          # Gapfill within the scaffold on the RHS of the new gap
          second_scaffold = scaffold.contigs[i+1]
          connections.push gapfill_scaffold.call(second_scaffold)
        end
      end
    end
  end

  def setup_output_directory(given_directory)
    output_directory = File.absolute_path(given_directory)
    log.debug "Using output directory: #{output_directory}" if log.debug?

    if File.exist?(output_directory)
      if !File.directory?(output_directory)
        log.error "Specified --output-directory #{output_directory} exists but is a file and not a directory. Cannot continue."
        exit 1
      elsif !File.writable?(output_directory)
        log.error "Specified --output-directory #{output_directory} is not writeable. Cannot continue."
        exit 1
      else
        log.debug "Already existing output directory #{output_directory} seems usable"
      end
    else
      # Creating a new output directory
      Dir.mkdir(output_directory)
    end

    return output_directory
  end

  def wander_a_genome(wanderer, genome, master_probed_graph, options)
    # Create new finishm_graph with only probes from the ends of the scaffolds of this genome
    probe_indices = []
    genome.each_numbered_probe{|probe| probe_indices.push(probe.number)}
    genome_graph = master_graph.subgraph(probe_indices)

    num_scaffolds = genome.scaffolds.length

    all_connections = wanderer.probed_graph_to_connections(genome_graph, scaffold_sequences, options)

    interpreter = Bio::FinishM::ConnectionInterpreter.new(all_connections, (0...num_scaffolds))
    connections = interpreter.doubly_single_contig_connections
    log.debug "Found #{connections.length} connections between contigs that can be used for scaffolding" if log.debug?

    return interpreter.scaffolds(connections)
  end

  class Genome
    attr_accessor :scaffolds, :filename, :numbered_probes

    def generate_numbered_probes(overhang, starting_probe_number)
      @numbered_probes = []
      @probe_number_to_scaffold_and_contig_and_side = {}

      current_probe_number = starting_probe_number
      overly_short_sequence_count = 0
      @scaffolds.each_with_index do |scaffold, scaffold_index|
        scaffold.contigs.each_with_index do |contig, contig_index|
          if contig.sequence.length < 2*overhang
            log.warn "Not attempting to make connections from overly short contig: it is the #{contig_index+1}th contig in scaffold `#{scaffold.name}' from the genome in `#{@filename}')"
            overly_short_sequence_count += 1
            nil
          else
            sequence = contig.sequence

            probe1 = NumberedProbe.new
            probe1.contig = contig
            probe1.number = current_probe_number; current_probe_number += 1
            probe1.side = :start
            fwd2 = Bio::Sequence::NA.new(sequence[0...options[:contig_end_length]])
            probe1.sequence = fwd2.reverse_complement.to_s

            probe2 = NumberedProbe.new
            probe2.contig = contig
            probe2.number = current_probe_number; current_probe_number += 1
            probe2.side = :end
            probe2.sequence = sequence[(sequence.length-options[:contig_end_length])...sequence.length]

            @numbered_probes[scaffold_index] ||= []
            @numbered_probes[scaffold_index][contig_index] = [probe1, probe2]

            @probe_number_to_scaffold_and_contig_and_side[probe1.number] = [scaffold, contig, :start]
            @probe_number_to_scaffold_and_contig_and_side[probe2.number] = [scaffold, contig, :end]
          end
        end
      end
      log.debug "Generated #{current_probe_number-starting_probe_number} probes for #{@filename}" if log.debug?
    end

    def number_of_probes
      @numbered_probes.flatten.length
    end

    def each_numbered_probe
      @numbered_probes.flatten.each do |probe|
        yield probe
      end
    end

    # Return true if probe number given is the probe at the beginning of the scaffold
    # or false if it is at the end. raise if unknown.
    def probe_at_start_of_scaffold?(probe_index)
      scaffold, contig, side = @probe_number_to_scaffold_and_contig_and_side[probe_index]
      if side == :start
        return true
      elsif side == :end
        return false
      else
        raise
      end
    end
  end

  class NumberedProbe
    attr_accessor :number, :contig, :side, :sequence
  end
end
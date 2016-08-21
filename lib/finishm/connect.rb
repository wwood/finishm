require 'tmpdir'

class Bio::FinishM::GapFiller
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm distance --read-sets <set1.fasta>,<set2.fasta>[,<set3.fasta>...] <assembly-specification>

Takes two or more set of reads determines the minimum distance between in a de-Bruijn graph traversal.

example: finishm distance --read-sets gene1reads.fa,gene2reads.fa,gene3reads.fa --fastq-gz reads.1.fq.gz,reads.2.fq.gz
\n"

    options.merge!({
      :graph_search_leash_length => 20000,
      })

    optparse_object.separator "\nRequired arguments:\n\n"
    optparse_object.on("--read-sets FILES", Array, "comma-separated list of fasta files containing read sets to connect [required]") do |arg|
      options[:read_sets] = arg
    end

    optparse_object.separator "\nThere must be some definition of of how to do the assembly, or else a path to a previous assembly directory:\n\n"
    Bio::FinishM::ReadInput.new.add_options(optparse_object, options)
    Bio::FinishM::GraphGenerator.new.add_options optparse_object, options

    optparse_object.separator "\nOptional graph search options:\n\n"
    optparse_object.on("--leash-length NUM", Integer, "Don't explore too far in the graph, only this many base pairs and not (much) more [default: #{options[:graph_search_leash_length] }]") do |arg|
      options[:graph_search_leash_length] = arg
    end
  end

  def validate_options(options, argv)
    #TODO: give a better description of the error that has occurred
    #TODO: require reads options
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      [
        :read_sets,
        ].each do |sym|
          if options[sym].nil?
            return "No option found to specify #{sym}"
          end
        end

      #if return nil from here, options all were parsed successfully
      return Bio::FinishM::ReadInput.new.validate_options(options, [])
    end
  end

  def run(options, argv)
    # read in fasta file of read sets

    # create a finishm graph with each of the reads in the read sets as probes

    # Determine which nodes contain the reads, and choose two reads from each detected node as examples

    # Dijkstra out from each of the probe nodes, ignoring direction

    # Detemine connections between readsets given the minimum distances

    # Output connectivity information
  end
end

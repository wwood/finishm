# A class representing reads or sets of reads to be assembled
class Bio::FinishM::ReadInput
  READ_INPUT_SYMBOLS = [
    :fasta_singles, :fastq_singles, :fasta_singles_gz, :fastq_singles_gz,
    #:interleaved_fasta, :interleaved_fastq, :interleaved_fasta_gz, :interleaved_fastq_gz,
    ]
  READ_INPUT_SYMBOLS.each do |sym|
    attr_accessor sym
  end

  # Given an OptionParser, add options to it, which parse out read-related options
  def add_options(option_parser, options)
    {
      '--fasta' => :fasta_singles,
      '--fastq' => :fastq_singles,
      '--fasta-gz' => :fasta_singles_gz,
      '--fastq-gz' => :fastq_singles_gz,
#       '--interleaved-fasta' => :interleaved_fasta,
#       '--interleaved-fastq' => :interleaved_fastq,
#       '--interleaved-fasta-gz' => :interleaved_fasta_gz,
#       '--interleaved-fastq-gz' => :interleaved_fastq_gz,
      }.each do |flag, sym|
        option_parser.on("#{flag} PATH", Array, "One or more paths to reads, comma separated") do |arg|
          options[sym] = arg
        end
      end
  end

  # Require at least 1 set of reads to be given, of any type
  def validate_options(options, argv)
    return nil if options[:previous_assembly] #bit of a hack, but hey
    READ_INPUT_SYMBOLS.each do |sym|
      return nil if options[sym]
    end
    return "No definition of reads for assembly was found"
  end

  # Parse options from options hash into instance variables for this object
  def parse_options(options)
    READ_INPUT_SYMBOLS.each do |sym|
      send("#{sym}=",options[sym]) if options[sym]
    end
  end

  # Output a string to be used on the command line with velvet
  def velvet_read_arguments
    readset_index = 1
    args = ''
    #Put paired sequences first in the hash (in Ruby, hashes are ordered) so that if they are paired, then odd numbered sequences
    # are read1 and even numbered sequences are read2
    {
#       :interleaved_fasta => '-fasta -shortPaired',
#       :interleaved_fastq => '-fastq -shortPaired',
#       :interleaved_fasta_gz => '-fasta.gz -shortPaired',
#       :interleaved_fastq_gz => '-fastq.gz -shortPaired',
      :fasta_singles => '-fasta -short',
      :fastq_singles => '-fastq -short',
      :fasta_singles_gz => '-fasta.gz -short',
      :fastq_singles_gz => '-fastq.gz -short',
      }.each do |sym, velvet_flag|
        paths = send(sym)
        unless paths.nil? or paths.empty?
          args += " #{velvet_flag}"
          paths.each do |path|
            args += " #{path}"
          end
        end
      end
    return args
  end
end

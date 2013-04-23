#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'
require 'progressbar'
require 'bio-ipcress'
require 'tempfile'

$:.unshift File.join(ENV['HOME'],'git','bioruby-primer3','lib')
require 'bio-primer3'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :logger_trace_level => 'info',
    :melting_temperature_optimum => nil,
    :melting_temperature_tolerance => 2,
    :min_primer_size => 15,
    :extra_global_primer3_options => {
      'PRIMER_MAX_POLY_X' => 4,
    },
    :persevere => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} [options]
      
      Takes a collection of contigs that are assumed to be a single circular genome.
      
      Design primers off the ends of them so that one big PCR might work.
      \n\n"
      
    opts.separator "Required arguments:\n\n"
    opts.on("-c", "--contigs FASTA_FILE", "A fasta file of contigs to be worked with [required]") do |arg|
      options[:contigs_file] = arg
    end
      
    opts.on("--min-distance-from-contig-ends DISTANCE", "Primers must be at least this far from the ends of the contigs [required]") do |arg|
      options[:min_distance] = arg.to_i
      raise Exception, "--min-distance-from-contig-ends has to be greater than/equal to 0, found #{arg}" unless options[:min_distance] >= 0
    end      
    opts.on("--max-distance-from-contig-ends DISTANCE", "Primers must be at most this far from the ends of the contigs [required]") do |arg|
      options[:max_distance] = arg.to_i
      raise Exception, "--max-distance-from-contig-ends has to be greater than 0, found #{arg}" unless options[:max_distance] > 0
    end
    
    opts.separator "\nOptional arguments:\n\n"
    opts.on("--optimum-melting-temperature TEMPERATURE", "Primers aim for this melting temperature [default: default in primer3 (currently 60C)]") do |arg|
      options[:melting_temperature_optimum] = arg.to_i
      raise Exception, " has to be greater than 0, found #{arg}" unless options[:melting_temperature_optimum] > 0
    end
    opts.on("--contig-universe FASTA_FILE", "All contigs in the mixture [default: unspecified (don't test this)]") do |arg|
      options[:contig_universe] = arg
    end
    
    opts.on("--persevere", "Don't automatically exit when a primer pair doesn't validate, though continue warning on ERROR log level [default: #{options[:persevere]}]") do
      options[:persevere] = true
    end
    
    opts.on("--primer3-options OPTION_LIST", "Give extra instructions to Primer3 [default <none>]. Acceptable values can be found in the primer3 manual e.g. 'PRIMER_MAX_POLY_X=4;PRIMER_MAX_SIZE=22' will specify those 2 parameters to primer3. Argument names are auto-capitalised so 'primer_max_poly_X=4;primer_max_size=22'is equivalent.") do |arg|
      options[:extra_global_primer3_options] = {}
      arg.split(';').each do |a2|
        splits = a2.split('=')
        unless splits.length == 2
          raise "Unexpected format of the --primer3-options flag, specifically couldn't parse this part: '#{a2}'"
        end
        options[:extra_global_primer3_options][splits[0].upcase]=splits[1]
      end
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:logger_trace_level] = 'error'}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:logger_trace_level] = s}
  end; o.parse!
  if ARGV.length != 0 or options[:min_distance].nil? or options[:max_distance].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); 
  Bio::Log::CLI.trace(options[:logger_trace_level])
  Bio::Log::CLI.configure('bio-primer3')
  Bio::Log::CLI.configure(LOG_NAME)
  # Setup logging for the primer3 gem too

  
  # Read the contigs in
  contigs = []
  Bio::FlatFile.foreach(options[:contigs_file]) do |entry|
    contigs.push entry
  end
  log.info "Read in #{contigs.length} contigs from #{options[:contigs_file]}"
  
  min_length = contigs.collect{|contig| contig.seq.length}.min
  log.info "Minimum contig length #{min_length}"
  unless options[:min_distance] < min_length/2
    log.error "Minimum primer distance from the ends of the contigs is too small, as the smallest contig is #{min_length} long, and the min distance must be at least twice this distance"
    exit 1
  end
  
  class PrimerList
    START_OF_CONTIG = 'start_of_contig'
    END_OF_CONTIG = 'end_of_contig'
    
    attr_accessor :contig_side, :primers, :contig_name
  end
  
  extra_primer3_options = {}
  unless options[:min_primer_size].nil?
    extra_primer3_options.merge!({
     'PRIMER_MIN_SIZE' => options[:min_primer_size],
    })
  end
  unless options[:melting_temperature_optimum].nil?
    extra_primer3_options.merge!({
      'PRIMER_OPT_TM' => options[:melting_temperature_optimum],
      'PRIMER_MIN_TM' => options[:melting_temperature_optimum]-options[:melting_temperature_tolerance],
      'PRIMER_MAX_TM' => options[:melting_temperature_optimum]+options[:melting_temperature_tolerance],
    })
  end
  unless options[:extra_global_primer3_options].nil?
    extra_primer3_options.merge! options[:extra_global_primer3_options]
  end
  if log.debug?
    # Get "debug-mode" from primer3 as well.
    extra_primer3_options.merge! 'PRIMER_EXPLAIN_FLAG' => '1'
  end
  
  # Predict a bunch of different primers for each end of each contig. Predict the start and end of each contig as the pair to pass to primer3
  primer3_results = []
  contigs.each do |contig|
    start_chunk = contig.seq[options[:min_distance]..options[:max_distance]]
    end_chunk = contig.seq[(contig.length-options[:max_distance]) .. (contig.length-options[:min_distance])].downcase
    log.debug "Start chunk length #{start_chunk.length}, end chunk length #{end_chunk.length}"
    
    # Join them together so that a forward primer will point off the end of the contig,
    # and a reverse primer will point off the start of the contig
    num_ns = 100
    joined = end_chunk+'N'*num_ns+start_chunk
    
    # Predict with primer3
    result = Bio::Primer3.run({
      'SEQUENCE_TEMPLATE' => joined,
      'PRIMER_TASK' => 'pick_sequencing_primers',
      'SEQUENCE_TARGET' => "#{end_chunk.length},#{num_ns}",
      'PRIMER_NUM_RETURN'=>'5',
      'PRIMER_PRODUCT_SIZE_RANGE'=>"#{num_ns}-#{joined.length}",
    }.merge(extra_primer3_options)
    )
    unless result.yeh?
      log.warn "No primers found for contig #{contig.definition}, giving up"
      exit 1
    end
    
    # Push each of the reported primers
    fwds = []
    reverses = []
    (0...result['PRIMER_LEFT_NUM_RETURNED'].to_i).each do |pair_number|
      fwds.push result["PRIMER_RIGHT_#{pair_number}_SEQUENCE"]
      reverses.push result["PRIMER_LEFT_#{pair_number}_SEQUENCE"]
    end
    contig_name = contig.definition
    
    f = PrimerList.new
    f.contig_side = PrimerList::START_OF_CONTIG
    f.primers = fwds
    f.contig_name = contig_name
    primer3_results.push f
    
    r = PrimerList.new
    r.contig_side = PrimerList::END_OF_CONTIG
    r.primers = reverses
    r.contig_name = contig_name
    primer3_results.push r
  end
  log.info "Finished getting first round of primers, now have #{primer3_results.length} sets of primers e.g. #{primer3_results[0].inspect}"
  
  # Sort from minimum number to most number of primers found
  # primer3_results.sort! do |a,b|
    # a.primers.length <=> b.primers.length
  # end
  # min_set = primer3_results[0]
  # log.info "Minimum number of primers found #{min_set.primers.length}, from #{min_set.contig_name}"

  if log.debug?
    log.debug "Primer sets that went in:"
    primer3_results.each do |res|
      log.debug res.inspect
    end
  end
  
  
  
  
  # Greedily try to find a set of primers such that one primer is picked from each location,
  # and this total set of primers doesn't conflict in any way
  
  # while not finished getting through the entire set
  # Progress is measured as a lsit of indices. Once the indices length is greater than
  failed = false
  current_path = [0]
  next_index_to_change = 0
  primer_sets = primer3_results.collect{|s| s.primers}
  
  while !failed and next_index_to_change < primer3_results.length
    if next_index_to_change == 0
      # garaunteed to be ok since there is only 1 primer
      current_path[0]=0
    
    # If there is no more, then we've failed.
    if current_path[0]==primer_sets.length
      failed = true
    else
      next_index_to_change += 1
    end
    
    else
      # Change the next possible index
      current_path[next_index_to_change] ||= -1
      if current_path[next_index_to_change] >= primer_sets[next_index_to_change].length
        # No more possibilities are available from this primer set, so have to backtrack
        next_index_to_change -= 1
      else
        current_path[next_index_to_change] += 1
        
        # Test whether this new primer conflicts with any of the old primers
        primer_to_test = primer_sets[next_index_to_change][current_path[next_index_to_change]]
        previous_primer_list = []
        (0...next_index_to_change).each do |i|
          previous_primer_list.push primer_sets[i][current_path[i]]
        end
        failed_against_prev = false
        previous_primer_list.each_with_index do |prev, i|
          log.debug "Testing #{primer_to_test.inspect} against #{prev.inspect}"
          if Bio::Primer3.test_primer_compatibility(primer_to_test, prev, extra_primer3_options) == false
            log.debug "Incompatible, unfortunately"
            log.error "This route through the code has never been tested, so you'll have to take a look at the code to ensure there is no bugs. Exiting."
            exit 1
            failed_against_prev = true
            break
          else
            log.debug "Compatible, cool"
          end
        end
        
        if failed_against_prev
          log.debug "At least one primer was incompatible, trying again"
          # Do nothing, try to change the same index again
        else
          log.debug "All compatible, cool. Now moving to the next index"
          current_path[next_index_to_change+1] = nil
          next_index_to_change += 1
        end
      end
    end
  end
  
  if failed
    log.error "Sorry, no sets of primers satisfy the criteria"
  else
    
    # First just make sure that everything is ok here
    log.info "Double checking to make sure there is no incompatibilities between primer pairs"
    num_compared = 0
    (0...primer_sets.length).to_a.combination(2) do |array|
      primer1 = primer_sets[array[0]][current_path[array[0]]]
      primer2 = primer_sets[array[1]][current_path[array[1]]]
      result = Bio::Primer3.test_primer_compatibility(primer1, primer2, extra_primer3_options)
      num_compared += 1
      
      if result == false
        log.error "Programming error!! There was supposed to be an OK path, but that path wasn't OK in the validation (#{primer1} and #{primer2}) were the problem"
        exit 1 unless options[:persevere]
      end
    end
    log.info "Validated #{num_compared} different pairs of primers, they don't seem to conflict with each other at all, according to primer3's check primers thing"
    
    # Check using in-silico PCR that all is ok
    # First, running ipcress on the contigs not joined together shouldn't yield any products
    log.debug "in-silico PCR: making sure there are no spurious primer pairings within the contigs themselves"
    ipcress_options = {:min_distance => 1, :max_distance => 10000, :mismatches => 0}
    num_compared = 0
    (0...primer_sets.length).to_a.combination(2) do |array|
      primer1 = primer_sets[array[0]][current_path[array[0]]]
      primer2 = primer_sets[array[1]][current_path[array[1]]]
      
      primer_set = Bio::Ipcress::PrimerSet.new primer1, primer2
      result = Bio::Ipcress.run primer_set, options[:contigs_file], ipcress_options
      num_compared += 1
      log.debug "Ipcress output:"
      log.debug result.inspect
      
      unless result.length == 0
        set1 = primer3_results[array[0]]
        set2 = primer3_results[array[1]]
        
        log.error "Unanticipated products generated from primer pair #{set1.contig_name}/#{set1.contig_side}/#{primer1} and #{set2.contig_name}/#{set2.contig_side}/#{primer2}. Sorry, fail."
        exit 1 unless options[:persevere]
      end
    end
    log.info "Validated #{num_compared} different pairs of primers so that unanticipated products are not formed according to iPCRess, and there doesn't seem to be any of those. Yey."
    
    # For each pair of primers, join together the corresponding contigs and validate that a sequencing product would eventuate
    num_compared = 0
    contigs_hash = {}
    contigs.each do |contig|
      contigs_hash[contig.definition] = contig.seq
    end
    log.debug "in-silico PCR: making sure expected products are generated..."
    (0...primer_sets.length).to_a.combination(2) do |array|
      set1 = primer3_results[array[0]]
      set2 = primer3_results[array[1]]
      
      primer1 = set1.primers[current_path[array[0]]]
      primer2 = set2.primers[current_path[array[1]]]
      
      log.debug "Testing #{primer1} and #{primer2} (indicies #{array[0]} and #{array[1]}, contigs #{set1.contig_name}/#{set1.contig_side} and #{set2.contig_name}/#{set2.contig_side})"
      primer_set = Bio::Ipcress::PrimerSet.new primer1, primer2
      Tempfile.open('ipcress') do |tempfile|
        # Have to correctly orient the template sequences
        s = PrimerList::START_OF_CONTIG
        e = PrimerList::END_OF_CONTIG
        f1 = contigs_hash[set1.contig_name]
        r1 = ' '+Bio::Sequence::NA.new(contigs_hash[set1.contig_name]).reverse_complement.to_s.upcase
        f2 = contigs_hash[set2.contig_name]
        r2 = ' '+Bio::Sequence::NA.new(contigs_hash[set2.contig_name]).reverse_complement.to_s.upcase
        
        seqs_ordered = case [set1.contig_side, set2.contig_side]
        when [s,s] then log.debug('case s s'); [r1,f2]
        when [s,e] then log.debug('case s e'); [r1,r2]
        when [e,s] then log.debug('case e s'); [f1,f2]
        when [e,e] then log.debug('case e e'); [f1,r2]
        end
        
        tempfile.puts '>test'
        tempfile.puts seqs_ordered[0]
        tempfile.puts 'N'*100
        tempfile.puts seqs_ordered[1]
        tempfile.close
        
        log.debug "Testing with iPCRess #{primer1} and #{primer2} from #{set1.contig_name} and #{set2.contig_name}, respectively"
        log.debug "Input fasta file"
        log.debug `cat #{tempfile.path} >/tmp/ta` if log.debug?
        log.debug "Fasta file input stats:"
        log.debug `seqstat #{tempfile.path}` if log.debug?
        
        
        results = Bio::Ipcress.run primer_set, tempfile.path, ipcress_options
        num_compared += 1
        unless results.length == 1
          if results.length == 0
            log.error "Anticipated products not generated in the hypothetical scenario of #{primer1} and #{primer2}, from #{set1.contig_name} and #{set2.contig_name}, respectively"
            exit 1 unless options[:persevere]
          else
            log.error "Too many PCR products generated in the hypothetical scenario of #{primer1} and #{primer2}, from #{set1.contig_name} and #{set2.contig_name}, respectively"
            log.error "Specifically, these PCR products were generated:"
            results.each do |res|
              log.error res.inspect
            end
            exit 1 unless options[:persevere]
          end
        end
      end
      num_compared += 1
    end
    log.info "Validated #{num_compared} different pairs of primers so that the anticipated products are formed according to iPCRess, and all seem to be as expected. Yey."
    
    
    # If a universe of possible contigs was given, do any primer pairs match up?
    num_compared = 0
    if options[:contig_universe]
      log.info "in-silico PCR: testing the contig universe. This could probably be sped up by only making a single call to ipcress, but oh well."
      num_compared = 0
      (0...primer_sets.length).to_a.combination(2) do |array|
        primer1 = primer_sets[array[0]][current_path[array[0]]]
        primer2 = primer_sets[array[1]][current_path[array[1]]]
        
        primer_set = Bio::Ipcress::PrimerSet.new primer1, primer2
        results = Bio::Ipcress.run primer_set, options[:contig_universe], ipcress_options
        num_compared += 1
        
        if results.length > 0
          log.warn "Found #{results.length} matches between #{primer1} and #{primer2} in the contig universe, expected none."
          exit 1 unless options[:persevere]
        end
        print '.' if log.info?
      end
      puts if log.info?
    else
      log.info "Not checking to see if primers match any other contigs not targeted are picked up by these primers, because no universe was specified."
    end
    log.debug "Tested #{num_compared} pairs of primers to see if anything else was hit, warnings above if any did so"
    
    
    log.info "Hoorah! Found an ok set of primers"
    puts %w(
      Contig
      Side
      Index
      Primer
    ).join("\t")
    current_path.each_with_index do |primer_index, primer_set_index|
      break if primer_set_index >= primer3_results.length
      
      end_info = primer3_results[primer_set_index]
      primer = end_info.primers[primer_index]
      puts [
        end_info.contig_name,
        end_info.contig_side,
        primer_index,
        primer,
      ].join("\t")
    end
  end
  
  
  
  
  
  
  
  
  
  
  
  
end #end if running as a script










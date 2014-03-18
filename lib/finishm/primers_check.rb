require 'optparse'
require 'tempfile'

require 'bio-logger'
require 'bio'
require 'progressbar'
require 'bio-ipcress'
$:.unshift File.join(ENV['HOME'],'git','bioruby-primer3','lib')
require 'bio-primer3'

class Bio::FinishM::Primers::Checker
  include Bio::FinishM::Logging

  def add_options(optparse_object, options)
    optparse_object.banner = "\nUsage: finishm primers_check --primers <primers_file> [options]

Check that each pair of primers in the primers fil is compatible.
\n\n"

    options.merge!({
      :logger => 'stderr',
      :logger_trace_level => 'info',
      :melting_temperature_optimum => 56,
      :melting_temperature_tolerance => 2,
      :extra_global_primer3_options => {
        'PRIMER_MAX_POLY_X' => 4,
        'PRIMER_EXPLAIN_FLAG' => '1',
        },
      :persevere => false,
      })

    optparse_object.separator "Required arguments:\n\n"
    optparse_object.on("-p", "--primers PRIMERS_FILE", String, "A file of primers, newline separated [required]") do |arg|
      options[:primers_file] = arg
    end

    optparse_object.separator "\nOptional arguments:\n\n"
    optparse_object.on("--contig-universe FASTA_FILE", String, "All contigs in the mixture [default: unspecified (don't test this)]") do |arg|
      options[:contig_universe] = arg
    end
    optparse_object.on("--optimum-melting-temperature TEMPERATURE", Integer, "Primers aim for this melting temperature [default: default in primer3 (currently #{options[:melting_temperature_optimum]}C)]") do |arg|
      options[:melting_temperature_optimum] = arg.to_i
      raise Exception, " has to be greater than 0, found #{arg}" unless options[:melting_temperature_optimum] > 0
    end
    optparse_object.on("--primer3-options OPTION_LIST", "Give extra instructions to Primer3 [default <none>]. Acceptable values can be found in the primer3 manual e.g. 'PRIMER_MAX_POLY_X=4;PRIMER_MAX_SIZE=22' will specify those 2 parameters to primer3. Argument names are auto-capitalised so 'primer_max_poly_X=4;primer_max_size=22'is equivalent.") do |arg|
      options[:extra_global_primer3_options] = {}
      arg.split(';').each do |a2|
        splits = a2.split('=')
        unless splits.length == 2
          raise "Unexpected format of the --primer3-options flag, specifically couldn't parse this part: '#{a2}'"
        end
        options[:extra_global_primer3_options][splits[0].upcase]=splits[1]
      end
    end
  end

  def validate_options(options, argv)
    if argv.length != 0
      return "Dangling argument(s) found e.g. #{argv[0] }"
    else
      [
        :primers_file,
        ].each do |sym|
          if options[sym].nil?
            return "No option found to specify #{sym}."
          end
        end
      return nil
    end
  end

  def run(options, argv)
    Bio::Log::CLI.configure('bio-primer3')
    found_error = false

    # Read the primers in
    primers = []
    File.foreach(options[:primers_file]) do |line|
      line.strip!
      next if line.empty?
      unless line.match(/^[atgc]+$/i)
        raise "Malformed primer sequence '#{line}' found in the file - I'm just after one primer per line, that's all."
      end
      primers.push line
    end
    log.info "Read in #{primers.length} primers e.g. #{primers[0] }"
    raise "Need at least 2 primers!" unless primers.length >= 2


    primer3_options = options[:extra_global_primer3_options]
    primer3_options.merge!({
      'PRIMER_OPT_TM' => options[:melting_temperature_optimum],
      'PRIMER_MIN_TM' => options[:melting_temperature_optimum]-options[:melting_temperature_tolerance],
      'PRIMER_MAX_TM' => options[:melting_temperature_optimum]+options[:melting_temperature_tolerance],
      })

    # Check to make sure there is no incompatibilities primer vs primer:
    num_compared = check_primer_compatibilities(primers, primer3_options)
    if num_compared == false
      log.error "Found at least one incompatible primer set, giving up."
      exit 1
    end
    log.info "Validated #{num_compared} different pairs of primers, they don't seem to conflict with each other at all, according to primer3's check primers thing"


    # Check the contig universe
    if options[:contig_universe]
      log.info "in-silico PCR: testing the contig universe. This could probably be sped up by only making a single call to ipcress, but oh well."
      num_universe_hits = check_contig_universe(primers, options[:contig_universe])
      if num_universe_hits > 0
        log.warn "Found #{results.length} matches between #{primer1} and #{primer2} in the contig universe, expected none."
      else
        log.info "No sequence appears to be amplified with any two primers in the background set of contigs (this is a good thing)."
      end
    else
      log.info "Not checking to see if primers match any other contigs not targeted are picked up by these primers, because no universe of contigs was specified."
    end


    if found_error
      log.warn "At least one problem detected in the primer checking process"
    else
      log.warn "No primer issues detected"
    end
  end



  # Each pair of primers should be compatible
  def check_primer_compatibilities(primers, primer3_options)
    # First just make sure that everything is ok here
    log.info "Double checking to make sure there is no incompatibilities between primer pairs"
    num_compared = 0
    primers.to_a.combination(2) do |array|
      primer1 = array[0]
      primer2 = array[1]

      compatible, result = Bio::Primer3.test_primer_compatibility(primer1, primer2, primer3_options, :return_result => true)
      num_compared += 1
      log.debug "Primer3 returned: #{result.inspect}" if log.debug?

      if compatible == false
        log.warn "Found an incompatibility between #{primer1} and #{primer2}"
        log.warn "Primer3 output was: #{result.inspect}"
        return false
      end
    end
    return num_compared
  end




  # Does any primer pairs amplify anything in a large set of background (universe) contigs?
  # They shouldn't.
  def check_contig_universe(primers, contig_universe)
    num_hits = 0
    num_compared = 0

    primers.to_a.combination(2) do |array|
      primer1 = array[0]
      primer2 = array[1]
      num_compared += 1

      primer_set = Bio::Ipcress::PrimerSet.new primer1, primer2
      results = Bio::Ipcress.run primer_set, options[:contig_universe], ipcress_options

      if results.length > 0
        num_hits += 1
        log.warn "Found #{results.length} matches between #{primer1} and #{primer2} in the contig universe, expected none."
      end
    end

    log.debug "Tested #{num_compared} pairs of primers to see if anything else was hit, warnings above if any did so"
  end
end

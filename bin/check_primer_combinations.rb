#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

$:.unshift File.join(ENV['HOME'],'git','bioruby-primer3','lib')
require 'bio-primer3'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -p1 <primer1> -f2 <primer_list_file>
      
      Uses primer3's \"check primers\" to find whether primers match against each other\n\n"
      
    opts.on("--primer1 PRIMER", "Primer on one side [required]") do |arg|
      options[:primer1] = arg
    end
    opts.on("--primers2 PRIMER_FILE", "A list of primers in a file, newline separated [required]") do |arg|
      options[:primers2_file] = arg
    end

    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in input data
  primers1 = [options[:primer1]]
  primers2 = File.open(options[:primers2_file]).read.split("\n").collect{|c| c.strip}
  log.info "Read in #{primers1.length} left primers and #{primers2.length} right primers e.g. #{primers1[0]} and #{primers2[0]}"
  
  goods = 0
  bads = 0
  failed_to_run = 0
  primers1.each do |primer1|
    primers2.each do |primer2|
      begin
        result, obj = Bio::Primer3.test_primer_compatibility primer1, primer2, 'PRIMER_EXPLAIN_FLAG'=>1
        
        puts [
          primer1, primer2, result, obj['PRIMER_LEFT_EXPLAIN'], obj['PRIMER_RIGHT_EXPLAIN']
        ].join "\t"
        
        if result
          goods += 1
        else
          bads += 1
        end
        
      rescue Exception => e
        failed_to_run += 1
      end
    end
  end
  log.info "Found #{goods} OK primer pairs and #{bads} not OK primer pairs"
  log.warn "#{failed_to_run} weren't checked by Primer3 because it failed to run" if failed_to_run > 0
end #end if running as a script
#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
reqiure 'csv'

SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = {
  :logger => 'stderr',
  :log_level => 'info',
}
o = OptionParser.new do |opts|
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>

    Takes a list of PCR primers that were put in several lanes (not all primers in all lanes), and a list of bands that were found, and decipher which bands are the result of which primer pairs, as best as possible\n\n"

  opts.on("--bands-file", "tsv file, with the band names as the first column, and the lane numbers that they appear in as the second column (comma separated) [required]") do |arg|
    options[:bands_file] = arg
  end
  opts.on("--primers-file", "tsv file, with the lane names as the first column, and the set of primers numbers that are in each lane as the second column (comma separated) [required]") do |arg|
    options[:primers_file] = arg
  end

  # logger options
  opts.separator "\nVerbosity:\n\n"
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {options[:log_level] = 'error'}
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| options[:log_level] = s}
end; o.parse!
if ARGV.length != 0 or options[:bands_file].nil? or options[:primers_file].nil?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]); Bio::Log::CLI.trace(options[:log_level]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)


# Read in the bands
bands_to_lanes = {}
CSV.foreach(options[:bands_file], :col_sep => "\t") do |row|
  raise "Malformed bands file in this line: #{row.inspect}" unless row.length == 2

  band_name = row[0]
  raise "Two bands were labeled the same way, as #{band_name.inspect}" if bands_to_lanes.key?(band_name)

  lanes_of_this_band = row[1].split(/[,\s]/).collect{|c| c.strip}
  bands_to_lanes[band_name] = lanes_of_this_band
end

# Read in the primer sets
lanes_to_primers = {}
CSV.foreach(options[:primers_file], :col_sep => "\t") do |row|
  raise "Malformed primers file in this line: #{row.inspect}" unless row.length == 2

  lane_name = row[0]
  raise "Two lanes were labeled the same way, as #{lane_name.inspect}" if lanes_to_primers.key?(lane_name)

  primers_of_this_band = row[1].split(/[,\s]/).collect{|c| c.strip}
  lanes_to_primers[lane_name] = primers_of_this_band
end


# Go through each pairing of primers. Which primer sets explain each band?
all_primers = lanes_to_primers.keys.flatten.sort
lanes = lanes_to_primers.keys
bands = bands_to_lanes.keys

bands_to_explaining_primer_pairs = {}

bands.each do |band|
  all_primers.combination(2).each do |array|
    primer1 = array.sort[0]
    primer2 = array.sort[1]

    band_agreeswith_this_primer_pair = true
    lanes.each do |lane|
      band_is_in_this_lane = bands_to_lanes[band].include?(lane)
      primers_here = lanes_to_primers[lane]
      if band and (!primers_here.include?(primer1) or !primers_here.include?(primer2))
        band_agreeswith_this_primer_pair = false
      end
      if !band and (primers_here.include?(primer1) and primers_here.include?(primer2))
        band_agreeswith_this_primer_pair = false
      end
    end

    if band_agreeswith_this_primer_pair
      bands_to_explaining_primer_pairs[band] ||= []
      bands_to_explaining_primer_pairs[band].push array
    end
  end

  puts [
    band,
    bands_to_explaining_primer_pairs[band].collect{|a| "(#{a.join(',')})"}
  ].join("\t")
end

require 'tempfile'
require 'rspec'
require 'pp'
require 'systemu'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb relative to the directory containing this test scripts, and the name of this tes script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
base = File.join File.dirname(__FILE__),'..'
path_to_script = "rdmd -I#{base}/../BioD/ #{base}/bin/read_selection_by_kmer.d --quiet"

# Re-build at the start
status, stdout, stderr = systemu "rdmd --build-only -I#{base}/../BioD/ #{base}/bin/read_selection_by_kmer.d"
raise stderr unless stderr == ""

describe script_under_test do
  it 'should scripting test ok' do
    reads = %w(>whitelist_me ATGCCCC >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC)
    blacklist = %w(ATGG)

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --blacklist #{blacklist_file.path} --reads #{reads_file.path} --kmer-coverage-target 10"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCC).join("\n")+"\n")
  end

  it 'should work without blacklist' do
    reads = %w(>whitelist_me ATGCCCC >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC)
    blacklist = %w()

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --reads #{reads_file.path} --kmer-coverage-target 10"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCC >blacklist_me ATGCATGG).join("\n")+"\n")
  end

  it 'should default to a low kmer coverage' do
    reads = %w(>whitelist_me ATGCCCC >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC)
    blacklist = %w()

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --reads #{reads_file.path}"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCC).join("\n")+"\n")
  end

  it 'should handle multiple whitelisted components' do
    reads = %w(>whitelist_me ATGCCCC >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC AAAA)
    blacklist = %w(ATGG)

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --blacklist #{blacklist_file.path} --reads #{reads_file.path} --kmer-coverage-target 10"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCC >ignore_me AAAAAAAA).join("\n")+"\n")
  end

  it 'should handle multiple whitelisted components in the same sequence' do
    reads = %w(>whitelist_me ATGCCCCAAAA >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC AAAA)
    blacklist = %w(ATGG)

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --blacklist #{blacklist_file.path} --reads #{reads_file.path}"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCCAAAA).join("\n")+"\n")
  end

  it 'should handle multiple whitelisted components in the same sequence and one in revcomp' do
    reads = %w(>whitelist_me ATGCCCCTTTT >blacklist_me ATGCATGG >ignore_me AAAAAAAA)
    whitelist = %w(ATGC AAAA)
    blacklist = %w(ATGG)

    reads_file = Tempfile.new('reads'); reads_file.puts reads.join("\n"); reads_file.close
    whitelist_file = Tempfile.new('whitelist'); whitelist_file.puts whitelist.join("\n"); whitelist_file.close
    blacklist_file = Tempfile.new('blacklist'); blacklist_file.puts blacklist.join("\n"); blacklist_file.close
    status, stdout, stderr = systemu "#{path_to_script} --whitelist #{whitelist_file.path} --blacklist #{blacklist_file.path} --reads #{reads_file.path}"
    raise stderr unless stderr == ""
    status.exitstatus.should eq(0)
    stdout.should eq(%w(>whitelist_me ATGCCCCTTTT).join("\n")+"\n")
  end
end

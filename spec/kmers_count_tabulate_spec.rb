require 'rspec'
require 'pp'
require 'systemu'
require 'tempfile'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/_spec/,'')
def assert_equal(e,o); o.should eq(e); end
path_to_script = File.join(File.dirname(__FILE__),'..','bin',script_under_test)



describe script_under_test do
  it 'should single file test' do
    Tempfile.open('spec') do |temp1|
      temp1.puts 'AAA 1'
      temp1.puts 'AAT 2'
      temp1.close

      status, stdout, stderr = systemu "#{path_to_script} #{temp1.path}"
      stderr.should eq("")
      status.exitstatus.should eq(0)
      stdout.should eq(["\t#{File.basename temp1.path}",
      "AAA\t1",
      "AAT\t2"].join("\n")+"\n")
    end
  end

  it 'should two file test' do
    Tempfile.open('spec') do |temp1|
      temp1.puts 'AAA 1'
      temp1.puts 'AAT 2'
      temp1.close

      Tempfile.open('spec') do |temp2|
        temp2.puts 'AAA 1'
        temp2.puts 'ATA 3'
        temp2.close

        status, stdout, stderr = systemu "#{path_to_script} #{temp1.path} #{temp2.path}"
        stderr.should eq("")
        status.exitstatus.should eq(0)
        stdout.should eq(["\t#{File.basename temp1.path}\t#{File.basename temp2.path}",
        "AAA\t1\t1",
        "AAT\t2\t0",
        "ATA\t0\t3"].join("\n")+"\n")
      end
    end
  end
end

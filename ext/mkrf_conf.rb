#(c) Copyright 2011 Raoul Bonnal. All Rights Reserved. Modified by Ben Woodcroft, 2014

# create Rakefile for shared library compilation



path = File.expand_path(File.dirname(__FILE__))

path_external = File.join(path, "../lib/external")

version = File.open(File.join(path_external,"VERSION"),'r')
Version = version.read
version.close

File.open(File.join(path,"Rakefile"),"w") do |rakefile|
rakefile.write <<-RAKE
require 'rbconfig'
require 'fileutils'
include FileUtils::Verbose
require 'rake/clean'

path = File.expand_path(File.dirname(__FILE__))
path_external = File.join(File.dirname(__FILE__), "../lib/external")

task :compile do
  cd(File.join(File.dirname(__FILE__),'src')) do
    #sh "patch -p1 < ../bioruby.patch"
    case RbConfig::CONFIG['host_os']
      when /linux/

        # Create library with default install params
        $stdout.puts "Making velvet shared library with default parameters"
        sh "make MAXKMERLENGTH=255 finishm velveth velvetg"
        $stdout.puts "Finished make finishm"
        shared_location = 'obj/shared'
        cp(File.join(shared_location,"libfinishm.so.1.0"), path_external)

      when /darwin/
        raise NotImplementedError, "possibly will work, but bio-velvet_underground is not tested on OSX"
      when /mswin|mingw/ then raise NotImplementedError, "bio-velvet_underground library is not available for Windows platform"
    end #case
  end #cd
end

task :clean do
  cd(File.join(path,'src')) do
    sh "make clean"
  end
  rm(File.join(path_external,"libvelvet.so.1.0"))
end

task :default => [:compile]

RAKE

end

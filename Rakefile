# encoding: utf-8

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "finishm"
  gem.homepage = "http://github.com/wwood/finishm"
  gem.license = "GPL-3.0+"
  gem.summary = %Q{Genome improvement and finishing with or without further sequencing effort}
  gem.description = %Q{De-novo assemblies generally only provide draft genomes. FinishM is aimed at improving these draft assemblies.}
  gem.email = "donttrustben near gmail.com"
  gem.authors = ["Ben J. Woodcroft"]
  # dependencies defined in Gemfile

  gem.extensions = "ext/mkrf_conf.rb"

  # by default, velvet as a git submodule is not included when making the gem
  # but we need it to be.
  gem.files.include "ext/src/src/*"
  gem.files.include "ext/src/Makefile"
  gem.files.include "ext/src/License"
  gem.files.include "ext/src/third-party/**/*"
end
Jeweler::RubygemsDotOrgTasks.new

require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

RSpec::Core::RakeTask.new(:rcov) do |spec|
  spec.pattern = 'spec/**/*_spec.rb'
  spec.rcov = true
end

task :default => :spec

require 'yard'
YARD::Rake::YardocTask.new

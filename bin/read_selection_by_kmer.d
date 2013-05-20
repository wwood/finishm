#!/usr/bin/env rdmd

import std.stdio;
import std.string;
import std.conv;
import std.getopt;

void main(string[] args){
  string whitelistFile, blacklistFile = null;
  bool help = false;
  getopt(args,
    "whitelist",  &whitelistFile,
    "blacklist",  &blacklistFile,
    "help|h",         &help,
  );
  if(help){
    writeln("my helpa");
  }
  if(whitelistFile == null){stderr.writeln("Need to specify a newline-separeted list of whitelisted kmers as --whitelist <file>");}

  //read in a text file of kmers that we wish to find, the whitelist
  string[] whitelistKmers = new string[];
  int i = 0;
  foreach(line; File(whitelistFile).byLine()){
  writeln(line);
  writeln(i);
    whitelistKmers[i] = to!string(line);
    i += 1;
  }
  stderr.writefln("Read in %i kmers as a whitelist, e.g. %s",whitelistKmers.length,whitelistKmers[0]);

  //if blacklistFile is specified, read in a list of kmers that are blacklisted, otherwise make an empty array

  //Iterate through the fastq reads given.
  //If they contain one of the blacklist kmers, then go to next
  //How many of each whitelist kmers are found (including in the reverse complement)?
  //Remove from the whitelist any kmers where the maximal allowed number of kmers has already been reached
}


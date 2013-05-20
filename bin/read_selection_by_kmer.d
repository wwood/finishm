#!/usr/bin/env rdmd

import std.stdio;
import std.string;
import std.conv;
import std.getopt;
import std.file;
import std.array;
import bio.core.fasta;

void main(string[] args){
  string whitelistFile, blacklistFile, fastaFile = null;
  bool help = false;
  getopt(args,
    "whitelist",  &whitelistFile,
    "blacklist",  &blacklistFile,
    "reads",  &fastaFile,
    "help|h",         &help,
  );
  if(help){
    writeln("my helpa");
  }
  else if(whitelistFile is null){stderr.writeln("Error: Need to specify a newline-separeted list of whitelisted kmers as --whitelist <file>");}
  else if(fastaFile is null){stderr.writeln("Error: Need to specify a fasta file of reads to work with as --reads <fasta_file>");}
  else {

  //read in a text file of kmers that we wish to find, the whitelist
  string[] whitelist = split(cast(string) read(whitelistFile));
  stderr.writeln("Read in ",whitelist.length," kmers as a whitelist, e.g. ",whitelist.front);

  //if blacklistFile is specified, read in a list of kmers that are blacklisted, otherwise make an empty array
  string[] blacklist = [];
  if(blacklistFile != null){
    blacklist = split(cast(string) read(blacklistFile));
    stderr.writeln("Read in ",blacklist.length," blacklisted kmers e.g. ",blacklist[0]);
  } else {
    stderr.writeln("No blacklisted kmers specified");
  }

  //Iterate through the fastq reads given.
  auto fastas = fastaRecords(fastaFile);
  writeln(fastas);
  writeln(fastas.front);
  writeln(fastas.front);
  writeln(fastas.front.sequence)

  //nucleotideSequence(window, true) //reverse complement too
  foreach(seq; fastas.popFront)
  if(fasta)
  //If they contain one of the blacklist kmers, then skip
  //How many of each whitelist kmers are found (including in the reverse complement)?
  //Remove from the whitelist any kmers where the maximal allowed number of kmers has already been reached
}
}

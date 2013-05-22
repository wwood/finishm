#!/usr/bin/env rdmd

import std.stdio;
import std.string;
import std.conv;
import std.getopt;
import std.file;
import std.array;
import bio.core.fasta;
import bio.core.sequence;
import std.algorithm;

void main(string[] args){
  string whitelistFile, blacklistFile, fastaFile = null;
  bool help, verbose, quiet = false;
  int targetPerKmer = 1;
  getopt(args,
    "whitelist",  &whitelistFile,
    "blacklist",  &blacklistFile,
    "reads",  &fastaFile,
    "kmer-coverage-target", &targetPerKmer,
    "verbose", &verbose,
    "quiet", &quiet,
    "help|h",         &help,
  );
  if(help){
    writeln("my helpa");
  }
  else if(whitelistFile is null){stderr.writeln("Error: Need to specify a newline-separeted list of whitelisted kmers as --whitelist <file>");}
  else if(fastaFile is null){stderr.writeln("Error: Need to specify a fasta file of reads to work with as --reads <fasta_file>");}
  else {
    if(verbose) quiet=false;

    //read in a text file of kmers that we wish to find, the whitelist
    string[] whitelist = split(cast(string) read(whitelistFile));
    if (verbose)
      stderr.writeln("Read in ",whitelist.length," kmers as a whitelist, e.g. ",whitelist.front);

    //if blacklistFile is specified, read in a list of kmers that are blacklisted, otherwise make an empty array
    string[] blacklist = [];
    if(blacklistFile != null){
      blacklist = split(cast(string) read(blacklistFile));
      if(verbose)
        stderr.writeln("Read in ",blacklist.length," blacklisted kmers e.g. ",blacklist[0]);
    } else {
      if(verbose)
        stderr.writeln("No blacklisted kmers specified");
    }

    int[] whitelistCounts;
    whitelistCounts.length = whitelist.length;

    //Iterate through the fastq reads given.
    auto fastas = fastaRecords(fastaFile);
    bool all_accounted_for = false;
    foreach(seq; fastas){
      if (verbose)
        stderr.writeln("Inspecting ", seq);
      //If they contain one of the blacklist kmers, then skip
      string fwd = seq.sequence;
      string rev = to!string(nucleotideSequence(fwd, true));

      //How many of each whitelist kmers are found (including in the reverse complement)?
      bool whitelisted = false;
      foreach(i, white; whitelist){
        //don't look for whitelisted kmers
        if (whitelistCounts[i] < targetPerKmer){
          if (std.string.indexOf(fwd, white) != -1 || std.string.indexOf(rev, white) != -1){
            if(verbose)
              stderr.writeln("Conforms to whitelist entry #",i);
            whitelisted = true;
            whitelistCounts[i] += 1;

            //check that all whitelist counts have not been overflowed.
            if(!quiet)
              stderr.writeln("kmer index ",i," now accounted for");
            if(minCount(whitelistCounts)[0] >= targetPerKmer){
              if(verbose)
                stderr.writeln("All whitelisted kmers now accounted for");
              all_accounted_for = true; //all done, no more fasta entries required
            }
          }
        }
      }
      if(!whitelisted) break;

      //I'm sure there is a faster way to search for an array of strings within a particular string, but eh for now.
      bool blacklisted = false;
      foreach(black; blacklist){
        if (std.string.indexOf(fwd, black) != -1 || std.string.indexOf(rev, black) != -1){
          //blacklisted kmer found
          blacklisted = true;
          break;
        }
      }
      if(blacklisted){
        if(verbose)
          stderr.writeln(fwd," contains a blacklisted kmer, not including this one");
        continue;
      } else {
        if(verbose)
          stderr.writeln(fwd," not blacklisted");
      }

      //print this sequence, as it is whitelisted and not blacklisted
      writeln(">", seq.header);
      writeln(fwd);

      if(all_accounted_for) break;
    }

    //output the number of kmers that were sufficiently covered
    int num_counted = count!("a >= b")(whitelistCounts, targetPerKmer);
    int num_not_counted = whitelistCounts.length - num_counted;
    if(!quiet)
      stderr.writeln("Found ",num_counted," from the whitelist as expected and ",num_not_counted," not enough times");
  }
}

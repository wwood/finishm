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
import std.container;
import std.c.stdlib;

void main(string[] args){
  string whitelistFile, blacklistFile, fastaFile = null;
  bool help, verbose, quiet, debugging = false;
  int targetPerKmer = 1;
  getopt(args,
    "whitelist",  &whitelistFile,
    "blacklist",  &blacklistFile,
    "reads",  &fastaFile,
    "kmer-coverage-target", &targetPerKmer,
    "verbose", &verbose,
    "debugging", &debugging,
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
    auto whites = split(cast(string) read(whitelistFile));
    auto whitelist = redBlackTree(whites);
    if (verbose)
      stderr.writeln("Read in ",whitelist.length," kmers as a whitelist, e.g. ",whitelist.front);

    //if blacklistFile is specified, read in a list of kmers that are blacklisted, otherwise make an empty array
    string[] blacks;
    auto blacklist = redBlackTree!string();
    if(blacklistFile != null){
      blacks = split(cast(string) read(blacklistFile));
      blacklist = redBlackTree(blacks);
      if(verbose)
        stderr.writeln("Read in ",blacklist.length," blacklisted kmers e.g. ",blacklist.front);
    } else {
      if(verbose)
        stderr.writeln("No blacklisted kmers specified");
    }

    int[string] whitelistCounts;
    foreach(white; whites){
      whitelistCounts[white] = 0;
    }

    // Find the minimum length of kmer being searched for
    int whitelistMinLength = map!"a.length"(whites).reduce!"a<b ? a : b";
    int whitelistMaxLength = map!"a.length"(whites).reduce!"a<b ? b : a";
    if (whitelistMinLength != whitelistMaxLength){
      stderr.writeln("Kmers must be of uniform length, but these ones weren't..");
      exit(1);
    }

    int blacklistMinLength = -1;
    if (verbose)
      stderr.writeln("Minimum length of kmer in whitelist is ",whitelistMinLength);
    if (blacklistFile){
      blacklistMinLength = map!"a.length"(blacks).reduce!"a<b ? a : b";
      if (verbose)
        stderr.writeln("Minimum length of kmer in whitelist is ",blacklistMinLength);
    }

    int num_reads_whitelisted = 0;
    int num_reads_blacklisted = 0;

    //Iterate through the fastq reads given.
    auto fastas = fastaRecords(fastaFile);
    bool all_accounted_for = false;
    int range_end;
    foreach(seq; fastas){
      if (verbose)
        stderr.writeln("Inspecting ", seq);
      //If they contain one of the blacklist kmers, then skip
      string fwd = seq.sequence;
      string rev = to!string(nucleotideSequence(fwd, true));

      range_end = fwd.length - whitelistMinLength;

      //How many of each whitelist kmers are found (including in the reverse complement)?
      bool whitelisted = false;
      foreach(i; 0 .. range_end){
        auto kmer = fwd[i .. (i+whitelistMinLength)];
        if (debugging)
          stdout.writeln("Whitelist inspecting kmer ",kmer);
        if (kmer in whitelist && whitelistCounts[kmer] < targetPerKmer){
          whitelisted = true;
          whitelistCounts[kmer] += 1;
          if (whitelistCounts[kmer] >= targetPerKmer){
            if(verbose)
              stderr.writeln("kmer index ",i," now accounted for");
            whitelist.removeKey(kmer);
            if (whitelist.empty()){
              if(verbose)
                stderr.writeln("All whitelisted kmers now accounted for");
              all_accounted_for = true; //all done, no more fasta entries required
            }
          }
        }
      }
      if(!whitelisted) continue;

      //I'm sure there is a faster way to search for an array of strings within a particular string, but eh for now.
      bool blacklisted = false;
      if (blacklistFile != null){
        foreach(i; 0 .. range_end){
          auto kmer = fwd[i .. (i+whitelistMinLength)];
          if (kmer in blacklist){
            //blacklisted kmer found
            blacklisted = true;
            break;
          }
        }
      }

      if(blacklisted){
        num_reads_blacklisted += 1;
        if(verbose)
          stderr.writeln(fwd," contains a blacklisted kmer, not including this one");
        continue;
      } else {
        if(verbose)
          stderr.writeln(fwd," not blacklisted");
      }

      //print this sequence, as it is whitelisted and not blacklisted
      num_reads_whitelisted += 1;
      writeln(">", seq.header);
      writeln(fwd);

      if(all_accounted_for) break;
    }

    //output the number of kmers that were sufficiently covered

    ulong num_counted = count!("a >= b")(whitelistCounts.values, targetPerKmer);
    ulong num_not_counted = whitelistCounts.length - num_counted;
    if(!quiet){
      stderr.writeln("Found ",num_counted," from the whitelist as expected and ",num_not_counted," not enough times");
      stderr.writeln("There were ",num_reads_whitelisted," reads output, and ",num_reads_blacklisted," reads blacklisted");
    }
  }
}

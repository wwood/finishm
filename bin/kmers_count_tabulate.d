#!/usr/bin/env rdmd

import std.stdio;
import std.conv;
import std.string;
import std.regex;
import std.getopt;

void main(string[] args){
  bool usePercentages = false;
  int minCount = 1;
  string trace = "info";
  getopt(args,
    "percentage", &usePercentages,
    "min-count",  &minCount,
    "trace",      &trace
  );

  //Create an array of open file handles, one for each argument given
  auto filenames = args[1 .. $];
  int[] totalCounts = new int[filenames.length];
  int kmerLength = 0;

  foreach(i, file; filenames){
    int count = 0;
    auto f = File(file);
    char[] line;
    while (f.readln(line)){
    //foreach(line; f.byLine()){
      if (kmerLength==0){
        kmerLength = 101;//to!int(indexOf(line, " "));
        stderr.writeln("Detected kmer length of ",kmerLength);
      }
      int thisCount = to!int(line[kmerLength+1 .. $-1]);
      //int thisCount = to!int(line[kmerLength+1 .. $]);
      count += thisCount;
      //if (count > 4000002){break;}
    }
    totalCounts[i] = count;
  }
  stderr.writeln("Finished ccounting total kmers, totals were: ",totalCounts);

  if(false){
  bool allFinished = false;
  bool[] finished = new bool[filenames.length];
  foreach (f; finished){f=false;}

  File[] files = new File[filenames.length];
  foreach(i; 0 .. files.length){
    files[i] = File(filenames[i]);
  }

  struct KmerCount {
    char[101] kmer;
    int count;
  }
  KmerCount[] currentRows = new KmerCount[files.length];
  char[] lineBuffer;
  foreach (i; 0..currentRows.length){
    //read in the line to teh buffer, so memory is not reallocated
    files[i].readln(lineBuffer);
    currentRows[i].kmer = chomp(files[i].readln(currentRows[i].kmer))[0..kmerLength];
    currentRows[i].count = to!int(line[kmerLength+1..$]);
  }

  //write headers
  enum ctr = ctRegex!(".*/(.+)");
  foreach(f; filenames){
    write("\t",match(f, ctr).captures[1]);
  }
  writeln();

  string[] toPrint = new string[filenames.length+1];
  int kmersCounted = 0;
  while (!allFinished){
    kmersCounted += 1;
    if (kmersCounted % (1024*1024) == 0){stderr.writeln("Processed ",kmersCounted," kmers.");}
    //Find the lowest kmer
    string lowestKmer = null;
    foreach (kc; currentRows){
      if (lowestKmer == null || kc.kmer < lowestKmer){
        lowestKmer = kc.kmer;
      }
    }

    //Go through each file, printing the number of this kmer found
    int totalObservations = 0;
    toPrint[0] = lowestKmer;
    foreach (i, kc; currentRows){
      if (kc.kmer == lowestKmer){
        totalObservations += kc.count;
        if (usePercentages){
          toPrint[i+1] = to!string(to!float(kc.count)/totalCounts[i]);
        } else {
          toPrint[i+1] = to!string(kc.count);
        }

        // Read a new line in, check if this file is finished
        auto line = files[i].readln;
        if (line == null){
          finished[i] = true;
          allFinished = true; //guilty until proven innocent
          foreach(f; finished){
            if (!f){
              allFinished = false;
            }
          }
          currentRows[i].kmer = null;
          currentRows[i].count = -1;
        } else {
          //Regular line to be read in
          currentRows[i].kmer = line[0..kmerLength];
          currentRows[i].count = to!int(line[kmerLength+1..$-1]);
        }
      } else {
        toPrint[i+1] = "0";
      }
    }
    if (totalObservations >= minCount){
      writeln(join(toPrint, "\t"));
    }
  }
}}

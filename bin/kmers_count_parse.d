#!/usr/bin/env rdmd

import std.stdio;
import std.csv;
import std.typecons;
import std.getopt;
import std.algorithm;



void main(string[] args)
{
  bool usePercentages = false;
  int minCount = 1;
  string trace = "info";
  getopt(args,
    "percentage", &usePercentages,
    "min-count",  &minCount,
    "trace",      &trace
  );

  auto kmersFile = File(args[1]);

  int lineCount = 0;
  foreach (line; kmersFile.byLine()) {
    lineCount += 1;
    if (lineCount % (1024*1024) == 0){
      stderr.writeln("Parsed ",lineCount, " lines");
    }
    auto reader = csvReader!(Tuple!(string,
      int, int, int, int, int,
      int, int, int, int, int,
      int, int, int, int, int,
      int, int, int, int, int,
      int, int, int, int, int,
      int, int, int, int, int,
      int, int))(line, ' ');
    foreach (record; reader) {
      /*writeln(record[0]);
      writeln(record[1]);
      writeln(record[1..32]);*/
      int[] range =  [1,2,3,4];
      int sum = 0;
      foreach(T; record[1..32]){
        sum += T;
      }
      if (sum >= minCount){
        writeln(line);
      }
    }
  }
}

# FinishM

FinishM attempts two things:

1. Improve draft genomes]by considering the computational problem to be about finishing, not assembly in the traditional sense.
2. If the genome cannot be finished with the current sequence data, FinishM can design and interpret a novel multiplex PCR finishing + next-gen sequencing experiment that vastly reduces the number of PCR reactions required to finish a genome, and automates interpretation of the results.

## Information on the table
Metagenome and isolate assemblers generate contigs from reads, but still leave valuable information on the table. FinishM exploits this information to improve/finish a draft genome without any further laboratory-based work.

In isolate genome sequencing, most/all contigs will constitute the vast majority of the genome being sequenced, but this information is ignored. Unlike a traditional assembler FinishM does not attempt to directly extend contigs, but instead focuses on connecting already assembled contigs.

FinishM has several modes:
* Determine which contig ends are connected in the assembly graph. See `finishm wander`.
* FinishM 'gapfills' (replaces N characters) using a graph-theoretic approach that appears to outperform current gapfilling programs. See `finishm gapfill`.
* Sometimes a human is better able to interpret an assembly graph than a machine. FinishM creates human interpretable graph visualisations that let humans solve assembly problems. See `finishm visualise`.

## A new finishing technique

FinishM combines multiplex PCR reactions and next-generation sequencing to simulataneously scaffold across and sequence the missing areas of genomes of interest. If you are interested in this please get in touch directly, as it is thus far not stable enough for large scale public consumption (not yet!).

## Installation

First, you'll need Ruby (FinishM is tested on 2.1). Then to install:
```sh
gem install finishm
```

FinishM also has some external dependencies:
* Velvet
* GraphViz

There are further external dependencies specific to the finishing PCR:
* exonerate
* primer3
* The D programming language environment
* BioD
Basically, the finishing PCR software probably won't work (yet) anywhere except the author's machine. I'm working on it.

## Usage
After installation, a listing of the modes and their usage:
```sh
finishm
```

## Citation

A manuscript describing the tools described here is currently in preparation. However, FinishM reuses code from velvet, exonerate/iPCRess, and primer3, so these tools may be worth citing.

## Copyright

Copyright (c) 2012-2014 Ben J. Woodcroft. See LICENSE.txt for
further details.


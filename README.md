__WARNING__! FinishM is very alpha software and not ready for prime time. There are many unfinished parts of it, and many bugs. Please use with care, and don't judge the authors too harshly.

# FinishM

FinishM attempts to improve draft genomes by considering the computational problem to be about finishing, not assembly in the traditional sense.

## A finishing approach to assembly
Metagenome and isolate assemblers generate contigs from reads, but still leave valuable information on the table. FinishM exploits this information to improve/finish a draft genome without any further laboratory-based work.

In even a moderately successful assembly, resultant contigs constitute the vast majority of the genome being sequenced, but this fact is ignored by assemblers. Unlike a traditional assembler FinishM does not attempt to directly extend contigs, but instead focuses on connecting already assembled contigs.

FinishM has several modes:
* Attempt to improve a genome. See `finishm roundup`. This mode fulfills both the `wander` and `gapfill` modes.
* Determine which contig ends are connected in the assembly graph. See `finishm wander`.
* FinishM 'gapfills' (replaces N characters) using a graph-theoretic approach that appears to outperform current gapfilling programs. See `finishm gapfill`.
* Sometimes a human is better able to interpret an assembly graph than a machine. FinishM creates human interpretable graph visualisations that let humans solve assembly problems. See `finishm visualise`.
* Some other experimental _de-novo_ (non-finishing) metagenome assembly techniques are implemented in `finishm assemble`.

## Installation

First, you'll need Ruby (FinishM is tested on 2.1). Then to install:
```sh
gem install finishm
```

FinishM also has some external dependencies:
* clustalo (for `gapfilling`/`roundup`)
* GraphViz (for the `visualise` mode)

## Usage
After installation, a listing of the modes and their usage:
```sh
finishm
```

## Developing
To hack on finishm:
```
git clone https://github.com/wwood/finishm.git
cd finishm
bundle install
git submodule update --init
cd ext/src
git checkout -b finishm origin/finishm #possibly this step is not required for newer versions of git
make MAXKMERLENGTH=255 finishm velveth velvetg
cp obj/shared/libfinishm.so.1.0 ../../lib/external/
cd ../..
./bin/finishm -h
```

## Citation

A manuscript describing the tools described here is currently in preparation. However, FinishM reuses code from velvet, clustalo and BioRuby, so these tools may be worth citing.

## Copyright

Copyright (c) 2012-2015 Ben J. Woodcroft. See LICENSE.txt for
further details.


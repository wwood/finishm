# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-
# stub: finishm 0.0.8 ruby lib
# stub: ext/mkrf_conf.rb

Gem::Specification.new do |s|
  s.name = "finishm"
  s.version = "0.0.8"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib"]
  s.authors = ["Ben J. Woodcroft"]
  s.date = "2016-08-01"
  s.description = "De-novo assemblies generally only provide draft genomes. FinishM is aimed at improving these draft assemblies."
  s.email = "donttrustben near gmail.com"
  s.executables = ["finishm"]
  s.extensions = ["ext/mkrf_conf.rb"]
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.md"
  ]
  s.files = [
    ".document",
    ".gitmodules",
    ".rspec",
    "Gemfile",
    "LICENSE.txt",
    "README.md",
    "Rakefile",
    "VERSION",
    "bin/finishm",
    "ext/mkrf_conf.rb",
    "ext/src/Makefile",
    "ext/src/src/allocArray.c",
    "ext/src/src/allocArray.h",
    "ext/src/src/autoOpen.c",
    "ext/src/src/autoOpen.h",
    "ext/src/src/binarySequences.c",
    "ext/src/src/binarySequences.h",
    "ext/src/src/concatenatedGraph.c",
    "ext/src/src/concatenatedGraph.h",
    "ext/src/src/concatenatedPreGraph.c",
    "ext/src/src/concatenatedPreGraph.h",
    "ext/src/src/correctedGraph.c",
    "ext/src/src/correctedGraph.h",
    "ext/src/src/dfib.c",
    "ext/src/src/dfib.h",
    "ext/src/src/dfibHeap.c",
    "ext/src/src/dfibHeap.h",
    "ext/src/src/dfibpriv.h",
    "ext/src/src/fib.c",
    "ext/src/src/fib.h",
    "ext/src/src/fibHeap.c",
    "ext/src/src/fibHeap.h",
    "ext/src/src/fibpriv.h",
    "ext/src/src/globals.h",
    "ext/src/src/graph.c",
    "ext/src/src/graph.h",
    "ext/src/src/graphReConstruction.c",
    "ext/src/src/graphReConstruction.h",
    "ext/src/src/graphStats.c",
    "ext/src/src/graphStats.h",
    "ext/src/src/graphStructures.h",
    "ext/src/src/kmer.c",
    "ext/src/src/kmer.h",
    "ext/src/src/kmerOccurenceTable.c",
    "ext/src/src/kmerOccurenceTable.h",
    "ext/src/src/kseq.h",
    "ext/src/src/locallyCorrectedGraph.c",
    "ext/src/src/locallyCorrectedGraph.h",
    "ext/src/src/passageMarker.c",
    "ext/src/src/passageMarker.h",
    "ext/src/src/preGraph.c",
    "ext/src/src/preGraph.h",
    "ext/src/src/preGraphConstruction.c",
    "ext/src/src/preGraphConstruction.h",
    "ext/src/src/probe_node_finder.c",
    "ext/src/src/probe_node_finder.h",
    "ext/src/src/readCoherentGraph.c",
    "ext/src/src/readCoherentGraph.h",
    "ext/src/src/readSet.c",
    "ext/src/src/readSet.h",
    "ext/src/src/readToNode.c",
    "ext/src/src/readToNode.h",
    "ext/src/src/recycleBin.c",
    "ext/src/src/recycleBin.h",
    "ext/src/src/roadMap.c",
    "ext/src/src/roadMap.h",
    "ext/src/src/run.c",
    "ext/src/src/run.h",
    "ext/src/src/run2.c",
    "ext/src/src/runReadToNode.c",
    "ext/src/src/scaffold.c",
    "ext/src/src/scaffold.h",
    "ext/src/src/shortReadPairs.c",
    "ext/src/src/shortReadPairs.h",
    "ext/src/src/splay.c",
    "ext/src/src/splay.h",
    "ext/src/src/splayTable.c",
    "ext/src/src/splayTable.h",
    "ext/src/src/tightString.c",
    "ext/src/src/tightString.h",
    "ext/src/src/utility.c",
    "ext/src/src/utility.h",
    "ext/src/third-party/zlib-1.2.3/ChangeLog",
    "ext/src/third-party/zlib-1.2.3/FAQ",
    "ext/src/third-party/zlib-1.2.3/INDEX",
    "ext/src/third-party/zlib-1.2.3/Makefile",
    "ext/src/third-party/zlib-1.2.3/Makefile.in",
    "ext/src/third-party/zlib-1.2.3/README",
    "ext/src/third-party/zlib-1.2.3/adler32.c",
    "ext/src/third-party/zlib-1.2.3/adler32.o",
    "ext/src/third-party/zlib-1.2.3/algorithm.txt",
    "ext/src/third-party/zlib-1.2.3/amiga/Makefile.pup",
    "ext/src/third-party/zlib-1.2.3/amiga/Makefile.sas",
    "ext/src/third-party/zlib-1.2.3/as400/bndsrc",
    "ext/src/third-party/zlib-1.2.3/as400/compile.clp",
    "ext/src/third-party/zlib-1.2.3/as400/readme.txt",
    "ext/src/third-party/zlib-1.2.3/as400/zlib.inc",
    "ext/src/third-party/zlib-1.2.3/compress.c",
    "ext/src/third-party/zlib-1.2.3/compress.o",
    "ext/src/third-party/zlib-1.2.3/configure",
    "ext/src/third-party/zlib-1.2.3/contrib/README.contrib",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/buffer_demo.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/mtest.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/read.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/test.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib-streams.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib-streams.ads",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib-thin.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib-thin.ads",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib.adb",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib.ads",
    "ext/src/third-party/zlib-1.2.3/contrib/ada/zlib.gpr",
    "ext/src/third-party/zlib-1.2.3/contrib/asm586/README.586",
    "ext/src/third-party/zlib-1.2.3/contrib/asm586/match.S",
    "ext/src/third-party/zlib-1.2.3/contrib/asm686/README.686",
    "ext/src/third-party/zlib-1.2.3/contrib/asm686/match.S",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/Makefile",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/README",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/blast.c",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/blast.h",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/test.pk",
    "ext/src/third-party/zlib-1.2.3/contrib/blast/test.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/delphi/ZLib.pas",
    "ext/src/third-party/zlib-1.2.3/contrib/delphi/ZLibConst.pas",
    "ext/src/third-party/zlib-1.2.3/contrib/delphi/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/delphi/zlibd32.mak",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib.build",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib.chm",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib.sln",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/AssemblyInfo.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/ChecksumImpl.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/CircularBuffer.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/CodecBase.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/Deflater.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/DotZLib.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/DotZLib.csproj",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/GZipStream.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/Inflater.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/DotZLib/UnitTests.cs",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/LICENSE_1_0.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/dotzlib/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/README",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/infback9.c",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/infback9.h",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/inffix9.h",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/inflate9.h",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/inftree9.c",
    "ext/src/third-party/zlib-1.2.3/contrib/infback9/inftree9.h",
    "ext/src/third-party/zlib-1.2.3/contrib/inflate86/inffas86.c",
    "ext/src/third-party/zlib-1.2.3/contrib/inflate86/inffast.S",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream/test.cpp",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream/zfstream.cpp",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream/zfstream.h",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream2/zstream.h",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream2/zstream_test.cpp",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream3/README",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream3/TODO",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream3/test.cc",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream3/zfstream.cc",
    "ext/src/third-party/zlib-1.2.3/contrib/iostream3/zfstream.h",
    "ext/src/third-party/zlib-1.2.3/contrib/masm686/match.asm",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/bld_ml64.bat",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/gvmat64.asm",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/gvmat64.obj",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/inffas8664.c",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/inffasx64.asm",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/inffasx64.obj",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx64/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/bld_ml32.bat",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/gvmat32.asm",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/gvmat32.obj",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/gvmat32c.c",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/inffas32.asm",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/inffas32.obj",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/mkasm.bat",
    "ext/src/third-party/zlib-1.2.3/contrib/masmx86/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/ChangeLogUnzip",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/Makefile",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/crypt.h",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/ioapi.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/ioapi.h",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/iowin32.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/iowin32.h",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/miniunz.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/minizip.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/mztools.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/mztools.h",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/unzip.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/unzip.h",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/zip.c",
    "ext/src/third-party/zlib-1.2.3/contrib/minizip/zip.h",
    "ext/src/third-party/zlib-1.2.3/contrib/pascal/example.pas",
    "ext/src/third-party/zlib-1.2.3/contrib/pascal/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/pascal/zlibd32.mak",
    "ext/src/third-party/zlib-1.2.3/contrib/pascal/zlibpas.pas",
    "ext/src/third-party/zlib-1.2.3/contrib/puff/Makefile",
    "ext/src/third-party/zlib-1.2.3/contrib/puff/README",
    "ext/src/third-party/zlib-1.2.3/contrib/puff/puff.c",
    "ext/src/third-party/zlib-1.2.3/contrib/puff/puff.h",
    "ext/src/third-party/zlib-1.2.3/contrib/puff/zeros.raw",
    "ext/src/third-party/zlib-1.2.3/contrib/testzlib/testzlib.c",
    "ext/src/third-party/zlib-1.2.3/contrib/testzlib/testzlib.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/untgz/Makefile",
    "ext/src/third-party/zlib-1.2.3/contrib/untgz/Makefile.msc",
    "ext/src/third-party/zlib-1.2.3/contrib/untgz/untgz.c",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/readme.txt",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/miniunz.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/minizip.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/testzlib.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/zlib.rc",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/zlibstat.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/zlibvc.def",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/zlibvc.sln",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc7/zlibvc.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/miniunz.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/minizip.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/testzlib.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/testzlibdll.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/zlib.rc",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/zlibstat.vcproj",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/zlibvc.def",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/zlibvc.sln",
    "ext/src/third-party/zlib-1.2.3/contrib/vstudio/vc8/zlibvc.vcproj",
    "ext/src/third-party/zlib-1.2.3/crc32.c",
    "ext/src/third-party/zlib-1.2.3/crc32.h",
    "ext/src/third-party/zlib-1.2.3/crc32.o",
    "ext/src/third-party/zlib-1.2.3/deflate.c",
    "ext/src/third-party/zlib-1.2.3/deflate.h",
    "ext/src/third-party/zlib-1.2.3/deflate.o",
    "ext/src/third-party/zlib-1.2.3/example",
    "ext/src/third-party/zlib-1.2.3/example.c",
    "ext/src/third-party/zlib-1.2.3/examples/README.examples",
    "ext/src/third-party/zlib-1.2.3/examples/fitblk.c",
    "ext/src/third-party/zlib-1.2.3/examples/gun.c",
    "ext/src/third-party/zlib-1.2.3/examples/gzappend.c",
    "ext/src/third-party/zlib-1.2.3/examples/gzjoin.c",
    "ext/src/third-party/zlib-1.2.3/examples/gzlog.c",
    "ext/src/third-party/zlib-1.2.3/examples/gzlog.h",
    "ext/src/third-party/zlib-1.2.3/examples/zlib_how.html",
    "ext/src/third-party/zlib-1.2.3/examples/zpipe.c",
    "ext/src/third-party/zlib-1.2.3/examples/zran.c",
    "ext/src/third-party/zlib-1.2.3/gzio.c",
    "ext/src/third-party/zlib-1.2.3/gzio.o",
    "ext/src/third-party/zlib-1.2.3/infback.c",
    "ext/src/third-party/zlib-1.2.3/infback.o",
    "ext/src/third-party/zlib-1.2.3/inffast.c",
    "ext/src/third-party/zlib-1.2.3/inffast.h",
    "ext/src/third-party/zlib-1.2.3/inffast.o",
    "ext/src/third-party/zlib-1.2.3/inffixed.h",
    "ext/src/third-party/zlib-1.2.3/inflate.c",
    "ext/src/third-party/zlib-1.2.3/inflate.h",
    "ext/src/third-party/zlib-1.2.3/inflate.o",
    "ext/src/third-party/zlib-1.2.3/inftrees.c",
    "ext/src/third-party/zlib-1.2.3/inftrees.h",
    "ext/src/third-party/zlib-1.2.3/inftrees.o",
    "ext/src/third-party/zlib-1.2.3/libz.a",
    "ext/src/third-party/zlib-1.2.3/make_vms.com",
    "ext/src/third-party/zlib-1.2.3/minigzip",
    "ext/src/third-party/zlib-1.2.3/minigzip.c",
    "ext/src/third-party/zlib-1.2.3/msdos/Makefile.bor",
    "ext/src/third-party/zlib-1.2.3/msdos/Makefile.dj2",
    "ext/src/third-party/zlib-1.2.3/msdos/Makefile.emx",
    "ext/src/third-party/zlib-1.2.3/msdos/Makefile.msc",
    "ext/src/third-party/zlib-1.2.3/msdos/Makefile.tc",
    "ext/src/third-party/zlib-1.2.3/old/Makefile.riscos",
    "ext/src/third-party/zlib-1.2.3/old/README",
    "ext/src/third-party/zlib-1.2.3/old/descrip.mms",
    "ext/src/third-party/zlib-1.2.3/old/os2/Makefile.os2",
    "ext/src/third-party/zlib-1.2.3/old/os2/zlib.def",
    "ext/src/third-party/zlib-1.2.3/old/visual-basic.txt",
    "ext/src/third-party/zlib-1.2.3/old/zlib.html",
    "ext/src/third-party/zlib-1.2.3/projects/README.projects",
    "ext/src/third-party/zlib-1.2.3/projects/visualc6/README.txt",
    "ext/src/third-party/zlib-1.2.3/projects/visualc6/example.dsp",
    "ext/src/third-party/zlib-1.2.3/projects/visualc6/minigzip.dsp",
    "ext/src/third-party/zlib-1.2.3/projects/visualc6/zlib.dsp",
    "ext/src/third-party/zlib-1.2.3/projects/visualc6/zlib.dsw",
    "ext/src/third-party/zlib-1.2.3/qnx/package.qpg",
    "ext/src/third-party/zlib-1.2.3/trees.c",
    "ext/src/third-party/zlib-1.2.3/trees.h",
    "ext/src/third-party/zlib-1.2.3/trees.o",
    "ext/src/third-party/zlib-1.2.3/uncompr.c",
    "ext/src/third-party/zlib-1.2.3/uncompr.o",
    "ext/src/third-party/zlib-1.2.3/win32/DLL_FAQ.txt",
    "ext/src/third-party/zlib-1.2.3/win32/Makefile.bor",
    "ext/src/third-party/zlib-1.2.3/win32/Makefile.emx",
    "ext/src/third-party/zlib-1.2.3/win32/Makefile.gcc",
    "ext/src/third-party/zlib-1.2.3/win32/Makefile.msc",
    "ext/src/third-party/zlib-1.2.3/win32/VisualC.txt",
    "ext/src/third-party/zlib-1.2.3/win32/zlib.def",
    "ext/src/third-party/zlib-1.2.3/win32/zlib1.rc",
    "ext/src/third-party/zlib-1.2.3/zconf.h",
    "ext/src/third-party/zlib-1.2.3/zconf.in.h",
    "ext/src/third-party/zlib-1.2.3/zlib.3",
    "ext/src/third-party/zlib-1.2.3/zlib.h",
    "ext/src/third-party/zlib-1.2.3/zutil.c",
    "ext/src/third-party/zlib-1.2.3/zutil.h",
    "ext/src/third-party/zlib-1.2.3/zutil.o",
    "finishm.gemspec",
    "lib/assembly/a_b_visualiser.rb",
    "lib/assembly/acyclic_connection_finder.rb",
    "lib/assembly/all_orfs.rb",
    "lib/assembly/bad_format_writer.rb",
    "lib/assembly/bam_probe_read_selector.rb",
    "lib/assembly/bubbly_assembler.rb",
    "lib/assembly/c_probe_node_finder.rb",
    "lib/assembly/connection_interpreter.rb",
    "lib/assembly/contig_printer.rb",
    "lib/assembly/coverage_based_graph_filter.rb",
    "lib/assembly/depth_first_search.rb",
    "lib/assembly/dijkstra.rb",
    "lib/assembly/fluffer.rb",
    "lib/assembly/graph_explorer.rb",
    "lib/assembly/graph_generator.rb",
    "lib/assembly/height_finder.rb",
    "lib/assembly/hybrid_velvet_graph.rb",
    "lib/assembly/input_genome.rb",
    "lib/assembly/kmer_coverage_based_path_filter.rb",
    "lib/assembly/node_finder.rb",
    "lib/assembly/oriented_node_trail.rb",
    "lib/assembly/paired_end_assembler.rb",
    "lib/assembly/paired_end_neighbour_finder.rb",
    "lib/assembly/probed_graph.rb",
    "lib/assembly/read_input.rb",
    "lib/assembly/read_to_node.rb",
    "lib/assembly/scaffold_breaker.rb",
    "lib/assembly/sequence_hasher.rb",
    "lib/assembly/single_coherent_paths_between_nodes.rb",
    "lib/assembly/single_coherent_wanderer.rb",
    "lib/assembly/single_ended_assembler.rb",
    "lib/assembly/velvet_c_binding.rb",
    "lib/assembly/velvet_graph_sequence_extractor.rb",
    "lib/external/VERSION",
    "lib/finishm/assemble.rb",
    "lib/finishm/explore.rb",
    "lib/finishm/finisher.rb",
    "lib/finishm/fluff.rb",
    "lib/finishm/gapfiller.rb",
    "lib/finishm/orfs_finder.rb",
    "lib/finishm/path_counter.rb",
    "lib/finishm/primers.rb",
    "lib/finishm/primers_check.rb",
    "lib/finishm/roundup.rb",
    "lib/finishm/sequence.rb",
    "lib/finishm/visualise.rb",
    "lib/finishm/wander.rb",
    "lib/kmer_abundance_pattern.rb",
    "lib/kmer_multi_abundance_file.rb",
    "lib/oligo_designer.rb",
    "lib/priner.rb",
    "spec/acyclic_connection_finder_spec.rb",
    "spec/all_orfs_spec.rb",
    "spec/assemble_spec.rb",
    "spec/bubbly_assembler_spec.rb",
    "spec/c_node_finder_spec.rb",
    "spec/connection_interpreter_spec.rb",
    "spec/contig_printer_spec.rb",
    "spec/coverage_based_graph_filter_spec.rb",
    "spec/data/6_3e4e5e6e.1vANME.bam",
    "spec/data/6_3e4e5e6e.1vANME.bam.bai",
    "spec/data/acyclic_connection_finder/1/probes.fa",
    "spec/data/acyclic_connection_finder/1/random1.fa",
    "spec/data/acyclic_connection_finder/1/random1.sammy.fa.gz",
    "spec/data/acyclic_connection_finder/1/random2.fa",
    "spec/data/acyclic_connection_finder/1/random2.sammy.fa.gz",
    "spec/data/assembly/1_simple_bubble_uneven_coverage/random3000.fa",
    "spec/data/assembly/1_simple_bubble_uneven_coverage/random3000.slightly_changed.fa",
    "spec/data/assembly/1_simple_bubble_uneven_coverage/reads_combined.fa.gz",
    "spec/data/assembly_visualiser/Contig_6_1_to_250.fa.kmers31",
    "spec/data/assembly_visualiser/Contig_7_1_to_250.fa.kmers31",
    "spec/data/assembly_visualiser/Graph",
    "spec/data/assembly_visualiser/start_kmers1",
    "spec/data/bands.csv",
    "spec/data/c_probe_node_finder/1/CnyUnifiedSeq",
    "spec/data/c_probe_node_finder/1/CnyUnifiedSeq.names",
    "spec/data/c_probe_node_finder/1/Graph2",
    "spec/data/c_probe_node_finder/1/LastGraph",
    "spec/data/c_probe_node_finder/1/Log",
    "spec/data/c_probe_node_finder/1/PreGraph",
    "spec/data/c_probe_node_finder/1/Roadmaps",
    "spec/data/c_probe_node_finder/1/contigs.fa",
    "spec/data/c_probe_node_finder/1/stats.txt",
    "spec/data/contig_printer/1/HOWTO_RECREATE",
    "spec/data/contig_printer/1/contigs.fa",
    "spec/data/contig_printer/1/seq.fa",
    "spec/data/contig_printer/1/seq.fa.svg",
    "spec/data/contig_printer/1/seq.fa.velvet/Graph2",
    "spec/data/contig_printer/1/seq.fa.velvet/LastGraph",
    "spec/data/contig_printer/1/seq.fa.velvet/Log",
    "spec/data/contig_printer/1/seq.fa.velvet/PreGraph",
    "spec/data/contig_printer/1/seq.fa.velvet/Roadmaps",
    "spec/data/contig_printer/1/seq.fa.velvet/Sequences",
    "spec/data/contig_printer/1/seq.fa.velvet/contigs.fa",
    "spec/data/contig_printer/1/seq.fa.velvet/stats.txt",
    "spec/data/contig_printer/1/seq.faVseq2_1to550.fa.bam",
    "spec/data/contig_printer/1/seq.faVseq2_1to550.fa.bam.bai",
    "spec/data/contig_printer/1/seq.node12.fa",
    "spec/data/contig_printer/1/seq1_1to550.fa",
    "spec/data/contig_printer/1/seq2_1to550.fa",
    "spec/data/contig_printer/1/seq2_1to550.fa.fai",
    "spec/data/explore/1/2seqs.sammy.fa",
    "spec/data/explore/1/HOWTO_RECREATE.txt",
    "spec/data/explore/1/a.fa",
    "spec/data/explore/1/seq1_and_a.fa",
    "spec/data/explore/1/seq2.fa",
    "spec/data/fluff/1/2seqs.sammy.fa",
    "spec/data/fluff/1/HOWTO_RECREATE.txt",
    "spec/data/fluff/1/seq1.fa",
    "spec/data/fluff/1/seq2.fa",
    "spec/data/gapfilling/1/reads.fa",
    "spec/data/gapfilling/1/trail_with_Ns.fa",
    "spec/data/gapfilling/1/velvetAssembly/Graph2",
    "spec/data/gapfilling/1/velvetAssembly/LastGraph",
    "spec/data/gapfilling/1/velvetAssembly/Log",
    "spec/data/gapfilling/1/velvetAssembly/PreGraph",
    "spec/data/gapfilling/1/velvetAssembly/Roadmaps",
    "spec/data/gapfilling/1/velvetAssembly/Sequences",
    "spec/data/gapfilling/1/velvetAssembly/contigs.fa",
    "spec/data/gapfilling/1/velvetAssembly/stats.txt",
    "spec/data/gapfilling/2/HOWTO_recreate",
    "spec/data/gapfilling/2/reference.fa",
    "spec/data/gapfilling/2/reference_part1.fa",
    "spec/data/gapfilling/2/reference_part2.fa",
    "spec/data/gapfilling/2/sammy_reads.fa.gz",
    "spec/data/gapfilling/2/with_gaps.fa",
    "spec/data/gapfilling/3/HOWTO_recreate",
    "spec/data/gapfilling/3/reads.fa.gz",
    "spec/data/gapfilling/3/reference_part1.fa",
    "spec/data/gapfilling/3/reference_part2.fa",
    "spec/data/gapfilling/3/with_gaps.fa",
    "spec/data/gapfilling/4/HOWTO_recreate",
    "spec/data/gapfilling/4/reads.fa.gz",
    "spec/data/gapfilling/5/HOWTO_RECREATE",
    "spec/data/gapfilling/5/answer.fna",
    "spec/data/gapfilling/5/gappy.fna",
    "spec/data/gapfilling/5/reads.fa",
    "spec/data/gapfilling/5/velvet51_3.5/LastGraph",
    "spec/data/gapfilling/5/velvet51_3.5/Sequences",
    "spec/data/gapfilling/6/random1.fa",
    "spec/data/gapfilling/6/random2.fa",
    "spec/data/gapfilling/6/random_sequence_length_2000",
    "spec/data/gapfilling/6/reads.random1.fa.gz",
    "spec/data/gapfilling/6/reads.random2.fa.gz",
    "spec/data/gapfilling/6/to_gapfill.fa",
    "spec/data/kmer_profile_to_assembly/multiple_abundance_file1.csv",
    "spec/data/kmers_count1.csv",
    "spec/data/kmers_count2.csv",
    "spec/data/out",
    "spec/data/positive_latching_pair.fa",
    "spec/data/primers.csv",
    "spec/data/read_selection_by_kmer/blacklist1.txt",
    "spec/data/read_selection_by_kmer/input.fasta",
    "spec/data/read_selection_by_kmer/whitelist1.txt",
    "spec/data/read_selection_by_kmer/whitelist2.txt",
    "spec/data/read_to_node/1_a_graph/HOWTO_RECREATE.txt",
    "spec/data/read_to_node/1_a_graph/LastGraph",
    "spec/data/read_to_node/1_a_graph/ReadToNode.bin",
    "spec/data/read_to_node/2_no_read256_or_259/HOWTO_RECREATE.txt",
    "spec/data/read_to_node/2_no_read256_or_259/LastGraph",
    "spec/data/read_to_node/2_no_read256_or_259/ReadToNode.bin",
    "spec/data/read_to_node/3_no_last_read/LastGraph",
    "spec/data/read_to_node/3_no_last_read/ReadToNode.bin",
    "spec/data/t/details.txt",
    "spec/data/t/details.txt.srt",
    "spec/data/t/location.txt",
    "spec/data/t/location.txt.srt",
    "spec/data/tweak/1_gap_then_unscaffolded/answer.fa",
    "spec/data/tweak/1_gap_then_unscaffolded/reads.fa.gz",
    "spec/data/tweak/1_gap_then_unscaffolded/scaffolds.fa",
    "spec/data/tweak/2_second_genome/answer2.fa",
    "spec/data/tweak/2_second_genome/reads.fa.gz",
    "spec/data/tweak/3_variant/answer.fa",
    "spec/data/tweak/3_variant/lesser_answer.fa",
    "spec/data/tweak/3_variant/reads.fa.gz",
    "spec/data/tweak/3_variant/with_gaps.fa",
    "spec/data/velvet_test_trails/Assem/Graph",
    "spec/data/velvet_test_trails/Assem/Graph2",
    "spec/data/velvet_test_trails/Assem/LastGraph",
    "spec/data/velvet_test_trails/Assem/Log",
    "spec/data/velvet_test_trails/Assem/PreGraph",
    "spec/data/velvet_test_trails/Assem/Roadmaps",
    "spec/data/velvet_test_trails/Assem/Sequences",
    "spec/data/velvet_test_trails/Assem/a.svg",
    "spec/data/velvet_test_trails/Assem/contigs.fa",
    "spec/data/velvet_test_trails/Assem/stats.txt",
    "spec/data/velvet_test_trails/node_fwds.fa",
    "spec/data/velvet_test_trails/node_seqs.fa",
    "spec/data/velvet_test_trails/nodes_fwd_rev.fa",
    "spec/data/velvet_test_trails/read1.fa",
    "spec/data/velvet_test_trails/reads.fa",
    "spec/data/velvet_test_trails_reverse/Assem/LastGraph",
    "spec/data/velvet_test_trails_reverse/Assem/a.svg",
    "spec/data/velvet_test_trails_reverse/reads_reversed.fa",
    "spec/data/visualise/1/LastGraph",
    "spec/data/visualise/2_paired_end/HOWTO_RECREATE.txt",
    "spec/data/visualise/2_paired_end/rand1.fa",
    "spec/data/visualise/2_paired_end/rand2.fa",
    "spec/data/visualise/2_paired_end/with_gaps.fa",
    "spec/data/visualise/2_paired_end/with_gaps.read_pairs.fa.gz",
    "spec/data/wander/1/random1.fa",
    "spec/data/wander/1/random1.sammy.fa",
    "spec/depth_first_search_spec.rb",
    "spec/dijkstra_spec.rb",
    "spec/explore_spec.rb",
    "spec/fluffer_spec.rb",
    "spec/gapfiller_spec.rb",
    "spec/graph_explorer_spec.rb",
    "spec/graph_generator_spec.rb",
    "spec/height_finder_spec.rb",
    "spec/kmer_abundance_pattern_spec.rb",
    "spec/kmer_coverage_based_path_filter_spec.rb",
    "spec/kmer_profile_finder_spec.rb",
    "spec/kmers_count_tabulate_spec.rb",
    "spec/oriented_node_trail_spec.rb",
    "spec/paired_end_neighbours_spec.rb",
    "spec/paths_between_nodes_spec.rb",
    "spec/priner_spec.rb",
    "spec/read_input_spec.rb",
    "spec/read_selection_by_kmer_spec.rb",
    "spec/read_to_node_spec.rb",
    "spec/roundup_spec.rb",
    "spec/scaffold_breaker_spec.rb",
    "spec/sequence_spec.rb",
    "spec/single_coherent_paths_between_nodes_spec.rb",
    "spec/single_coherent_wanderer_spec.rb",
    "spec/single_ended_assembler_spec.rb",
    "spec/spec_helper.rb",
    "spec/velvet_graph_sequence_extractor_spec.rb",
    "spec/visualise_spec.rb",
    "spec/wander_spec.rb",
    "spec/watch_for_changes.sh",
    "validation/fasta_compare.rb",
    "validation/gapfill_simulate_perfect.rb"
  ]
  s.homepage = "http://github.com/wwood/finishm"
  s.licenses = ["GPL3+"]
  s.rubygems_version = "2.5.1"
  s.summary = "Genome improvement and finishing with or without further sequencing effort"

  if s.respond_to? :specification_version then
    s.specification_version = 4

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio-ipcress>, ["~> 0.0"])
      s.add_runtime_dependency(%q<bio-logger>, ["~> 1.0"])
      s.add_runtime_dependency(%q<bio>, [">= 1.4.3", "~> 1.4"])
      s.add_runtime_dependency(%q<progressbar>, ["~> 0.21"])
      s.add_runtime_dependency(%q<ruby-graphviz>, [">= 1.0.9", "~> 1.0"])
      s.add_runtime_dependency(%q<ds>, [">= 0.0.4", "~> 0.0"])
      s.add_runtime_dependency(%q<hopcsv>, [">= 0.4.3", "~> 0.4"])
      s.add_runtime_dependency(%q<bio-velvet>, ["~> 0.6"])
      s.add_runtime_dependency(%q<bio-velvet_underground>, ["~> 0.3"])
      s.add_runtime_dependency(%q<ruby-progressbar>, [">= 1.4.2", "~> 1.4"])
      s.add_runtime_dependency(%q<yargraph>, [">= 0.0.4", "~> 0.0"])
      s.add_runtime_dependency(%q<pry>, ["~> 0.10"])
      s.add_development_dependency(%q<rspec>, [">= 2.8.0", "~> 2.8"])
      s.add_development_dependency(%q<yard>, ["~> 0.7"])
      s.add_development_dependency(%q<rdoc>, ["~> 3.12"])
      s.add_development_dependency(%q<bundler>, [">= 1.6.2", "~> 1.6"])
      s.add_development_dependency(%q<jeweler>, [">= 2.0.1", "~> 2.0"])
      s.add_development_dependency(%q<bio-commandeer>, ["~> 0.1"])
    else
      s.add_dependency(%q<bio-ipcress>, ["~> 0.0"])
      s.add_dependency(%q<bio-logger>, ["~> 1.0"])
      s.add_dependency(%q<bio>, [">= 1.4.3", "~> 1.4"])
      s.add_dependency(%q<progressbar>, ["~> 0.21"])
      s.add_dependency(%q<ruby-graphviz>, [">= 1.0.9", "~> 1.0"])
      s.add_dependency(%q<ds>, [">= 0.0.4", "~> 0.0"])
      s.add_dependency(%q<hopcsv>, [">= 0.4.3", "~> 0.4"])
      s.add_dependency(%q<bio-velvet>, ["~> 0.6"])
      s.add_dependency(%q<bio-velvet_underground>, ["~> 0.3"])
      s.add_dependency(%q<ruby-progressbar>, [">= 1.4.2", "~> 1.4"])
      s.add_dependency(%q<yargraph>, [">= 0.0.4", "~> 0.0"])
      s.add_dependency(%q<pry>, ["~> 0.10"])
      s.add_dependency(%q<rspec>, [">= 2.8.0", "~> 2.8"])
      s.add_dependency(%q<yard>, ["~> 0.7"])
      s.add_dependency(%q<rdoc>, ["~> 3.12"])
      s.add_dependency(%q<bundler>, [">= 1.6.2", "~> 1.6"])
      s.add_dependency(%q<jeweler>, [">= 2.0.1", "~> 2.0"])
      s.add_dependency(%q<bio-commandeer>, ["~> 0.1"])
    end
  else
    s.add_dependency(%q<bio-ipcress>, ["~> 0.0"])
    s.add_dependency(%q<bio-logger>, ["~> 1.0"])
    s.add_dependency(%q<bio>, [">= 1.4.3", "~> 1.4"])
    s.add_dependency(%q<progressbar>, ["~> 0.21"])
    s.add_dependency(%q<ruby-graphviz>, [">= 1.0.9", "~> 1.0"])
    s.add_dependency(%q<ds>, [">= 0.0.4", "~> 0.0"])
    s.add_dependency(%q<hopcsv>, [">= 0.4.3", "~> 0.4"])
    s.add_dependency(%q<bio-velvet>, ["~> 0.6"])
    s.add_dependency(%q<bio-velvet_underground>, ["~> 0.3"])
    s.add_dependency(%q<ruby-progressbar>, [">= 1.4.2", "~> 1.4"])
    s.add_dependency(%q<yargraph>, [">= 0.0.4", "~> 0.0"])
    s.add_dependency(%q<pry>, ["~> 0.10"])
    s.add_dependency(%q<rspec>, [">= 2.8.0", "~> 2.8"])
    s.add_dependency(%q<yard>, ["~> 0.7"])
    s.add_dependency(%q<rdoc>, ["~> 3.12"])
    s.add_dependency(%q<bundler>, [">= 1.6.2", "~> 1.6"])
    s.add_dependency(%q<jeweler>, [">= 2.0.1", "~> 2.0"])
    s.add_dependency(%q<bio-commandeer>, ["~> 0.1"])
  end
end


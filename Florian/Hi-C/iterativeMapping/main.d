// Iterative fragment mapping
// (c) Danny Arends (HU-Berlin) Sept - 2019

import std.stdio : File, writeln, writefln;
import std.typecons : refCounted;
import std.string : format, split, strip, toUpper, indexOf, splitLines;
import std.conv : to;
import std.zlib : Compress, compress, HeaderFormat, Z_SYNC_FLUSH, Z_FULL_FLUSH;
import std.file : exists, remove;
import std.process : executeShell;
import std.path : baseName;
import std.array : Appender;

// imports from external dub packages
import std.io : IOFile = File;
import containers.hashmap : HashMap;
import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip, zip;
import iopipe.bufpipe : bufd;

enum MB10 = 10 * 1024 * 1024;

struct MapInfo {
  size_t nAlignments;
  string alignLine;
}

struct AlignmentResult {
  size_t mapped;
  size_t unknown;
  string[] toremove;
}

struct FastQ{
  string path;
  size_t nlines;
  size_t nReads;
  size_t readLength;
  size_t nUnmapped;
  HashMap!(string, MapInfo) reads;
}

// Get basic information about the content of a gz fastq file
FastQ infoFastQ(string fastqpath, bool verbose = true) {
  FastQ fq = FastQ(fastqpath);
  auto fp = IOFile(fq.path).refCounted.bufd.unzip(CompressionFormat.gzip);
  string shortname;
  
  foreach (line; fp.assumeText.byLineRange) {
    final switch (fq.nlines % 4) {
      // Fasta header line
      case 0: fq.nReads++;
              fq.nUnmapped++;
              shortname = (to!string(line[1 .. $])).split(" ")[0];
              fq.reads[shortname.idup] = MapInfo(0, ""); // Read has not been aligned yet
              break;
      // Read
      case 1: fq.readLength = line.length; break;
      case 2: break; // +
      case 3: break; // Quality
    }
    fq.nlines++;
  }
  writefln("FastQ contained %s reads, read length: %s bp", fq.nReads, fq.readLength);
  return(fq);
}

// Reduce the readlength of the unmapped reads in a fastq file to the specified readLength
string reduceFastQ(ref FastQ fq, string fmt = "readlength%s.fq.gz", size_t readLength = 25) {
  size_t nlines = 0;
  string outputfilename = format(fmt, readLength);
  writefln("Writing reduced fastQ file to %s", outputfilename);
  
  auto fp = IOFile(fq.path).refCounted.bufd.unzip(CompressionFormat.gzip);
  auto ofp = File(outputfilename, "w");

  Compress cmp = new Compress(HeaderFormat.gzip);
  bool notAlignedYet = true;
  auto outputbuffer = Appender!(char[])();
  foreach (line; fp.assumeText.byLineRange!true) {
    if (nlines % 4 == 0) {
      string shortname = (to!string(line[1 .. $])).split(" ")[0];
      notAlignedYet = ((shortname in fq.reads) !is null);
    }
    if (nlines % 4 == 1 || nlines % 4 == 3) { // Read and Quality score need to be reduced
      if (notAlignedYet) {
        outputbuffer.put(line[0..readLength]);
        outputbuffer.put("\n");
      }
    } else {
      if (notAlignedYet) { outputbuffer.put(line); }
    }
    nlines++;
    if (outputbuffer.data.length > MB10) { // if we have over 10mb of buffer compress it and write it to the disk
      ofp.rawWrite(cmp.compress(outputbuffer.data));
      ofp.rawWrite(cmp.flush(Z_FULL_FLUSH)); // Flush to be sure
      outputbuffer.clear();
    }
  }
  ofp.rawWrite(cmp.compress(outputbuffer.data));
  ofp.rawWrite(cmp.flush());
  return(outputfilename);
}

// Map a fastq file to the reference genome and update the number of unmapped reads (nUnmapped)
void mapToGenome(ref FastQ fq, string fastqpath, string referencepath, string outputfilename = "alignments.txt", size_t minMapQ = 30) {
  // Align reads using BWA
  auto ofp = File(outputfilename, "a");
  string cmd = format("~/Github/bwa/bwa mem -v 2 -t 12 -T 10 %s %s", referencepath, fastqpath);
  writefln("Aligning reads from %s to %s using bwa", fastqpath, baseName(referencepath));
  auto ret = executeShell(cmd);
  writefln("output length: %s", ret.output.length);
  size_t lines = fq.processBWA(ret.output, minMapQ); // Process Output
  auto res = fq.parseAlignments(ofp); // Find mapped reads
  fq.removeFromAA(res.toremove); // Remove the mapped reads

  writefln("Parsed %s lines, mapped %s reads, unmapped reads left %s", lines, res.mapped, fq.nUnmapped);
}

// Process BWA output and update the fq.reads structure
size_t processBWA(ref FastQ fq, const string output, size_t minMapQ = 30){
  writefln("Processing %.2f Megabyte of BWA output", to!float(output.length) / (1024 * 1024));
  size_t lines = 0;
  size_t unknownreads = 0;
  foreach (line; output.splitLines()) {
    if (line[0] == '@') continue; // SAM header (just ignore it)
    lines++;
    auto sline = line.split("\t");
    if (sline.length > 4) { // Less than 4, not a valid SAM alignment line
      if (sline[2] == "*") continue; // Not aligned, just continue with the next line
      string qname = sline[0];
      if ((qname in fq.reads) is null) {
        unknownreads++;
        continue;
      }
      size_t mapq = to!size_t(sline[4]);
      // Good alignment (above minMapQ), if an alignment of the read was already found the mapQ doesn't matter
      if (fq.reads[qname].nAlignments > 0 || mapq > minMapQ) {
        fq.reads[qname] = MapInfo(fq.reads[qname].nAlignments + 1, line.idup);
      }
    }
  }
  writefln("Processed %d lines of BWA output (%s unknown reads)", lines, unknownreads);
  return(lines);
}

// Parse alignment results of BWA to figure out which keys need to be removed since they have 1 unique alignment
AlignmentResult parseAlignments(ref FastQ fq, File ofp) {
  AlignmentResult res;
  foreach (key; fq.reads.byKey) {
    if ((key in fq.reads) is null) {
      res.unknown++;
      continue;
    }
    if (fq.reads[key].nAlignments == 1) {
      ofp.write(fq.reads[key].alignLine ~ "\n");
      res.toremove ~= key;
      res.mapped++;
    } else { // none or multiple alignments found, reset the nAlignments to 0
      fq.reads[key] = MapInfo(0);
    }
  }
  return(res);
}

// Remove the keys from the array, return the number of keys not found in the array
size_t removeFromAA(ref FastQ fq, const string[] toremove) {
  size_t unknownreads = 0;
  foreach (key; toremove) {
    if ((key in fq.reads) is null) {
      unknownreads++;
      continue;
    }
    fq.reads.remove(key);
    fq.nUnmapped = fq.nUnmapped - 1;
  }
  return(unknownreads);
}

// dub -- /halde/Hi-C/Human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /halde/Hi-C/Human/ENCFF319AST.1Mio.fastq.gz ENCFF319AST.1Mio.alignment
// dub -- /halde/Hi-C/Human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /halde/Hi-C/Human/ENCFF478EAB.fastq.gz ENCFF478EAB.alignment
int main (string[] args) {
  if (args.length < 4) {
    writeln("Please provide the fastq input file and fasta reference");
    return(-1);
  }
  string referencepath = args[1]; //"/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa.gz";
  string fastqpath = args[2]; //"/halde/Hi-C/Human/ENCFF319AST.1Mio.fastq.gz";
  string outputfilename = args[3]; //"ReadAlignments.txt";
  string tmpfmt = baseName(args[3], ".alignment") ~ ".%s.fq.gz";
  writefln("tmpfmt: %s", tmpfmt);
  size_t readLength = 25;
  if (outputfilename.exists) {
    writefln("Deleting previous output file: %s", outputfilename);
    outputfilename.remove();
  }
  FastQ fq = infoFastQ(fastqpath);

  while (readLength <= fq.readLength) {
    if(fq.nUnmapped == 0) break;
    string path = fq.reduceFastQ(tmpfmt, readLength);
    fq.mapToGenome(path, referencepath, outputfilename);
    return(-1);
    if (path.exists) {
      writefln("Deleting temporary input file: %s", path);
      path.remove();
    }
    readLength += 5;
  }
  return(0);
}

// Iterative fragment mapping 

import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip, zip;
import iopipe.bufpipe : bufd;

import std.io : IOFile = File;
import std.stdio : File, writeln, writefln;
import std.typecons : refCounted;
import std.string : format, split, strip, toUpper, indexOf, splitLines;
import std.conv : to;
import std.zlib : Compress, compress, HeaderFormat;
import std.file : exists, remove;
import std.process : executeShell;
import std.path : baseName;

struct MapInfo {
  size_t nAlignments;
  string alignLine;
}

struct FastQ{
  string path;
  size_t nlines;
  size_t nReads;
  size_t readLength;
  size_t nUnmapped;
  MapInfo[string] reads;
}

bool has(const MapInfo[string] aa, string k){ return((k in aa)!is null); }

// Get basic information about the content of a gz fastq file
FastQ infoFastQ(string fastqpath) {
  FastQ fq = FastQ(fastqpath);
  auto fp = IOFile(fq.path).refCounted.bufd.unzip(CompressionFormat.gzip);

  foreach (line; fp.assumeText.byLineRange) {
    final switch (fq.nlines % 4) {
      // Fasta header line
      case 0: fq.nReads++;
              fq.nUnmapped++;
              string shortname = (to!string(line[1 .. $])).split(" ")[0];
              fq.reads[shortname] = MapInfo(0, ""); // Read has not been aligned yet
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
string reduceFastQ(const FastQ fq, size_t readLength = 25) {
  size_t nlines = 0;
  string outputfilename = format("readlength%s.fq.gz", readLength);
  writefln("Writing reduced fastQ file to %s", outputfilename);
  
  auto fp = IOFile(fq.path).refCounted.bufd.unzip(CompressionFormat.gzip);
  auto ofp = File(outputfilename, "w");
  Compress cmp = new Compress(HeaderFormat.gzip);
  bool notAlignedYet;
  foreach (line; fp.assumeText.byLineRange!true) {
    if (nlines % 4 == 0) {
      string shortname = (to!string(line[1 .. $])).split(" ")[0];
      if (fq.reads.has(shortname)) {
        notAlignedYet = true;
      } else {
        notAlignedYet = false;
      }
    }
    if (nlines % 4 == 1 || nlines % 4 == 3) { // Read and Quality score need to be reduced
      if (notAlignedYet) ofp.rawWrite(cmp.compress(line[0..readLength] ~ "\n"));
    } else {
      if (notAlignedYet) ofp.rawWrite(cmp.compress(line));
    }
    nlines++;
  }
  ofp.rawWrite(cmp.flush());
  return(outputfilename);
}

// Map a fastq file to the reference genome and update the number of unmapped reads (nUnmapped)
void mapToGenome(ref FastQ fq, string fastqpath, string referencepath, string outputfilename = "alignments.txt", size_t minMapQ = 30) {
  auto ofp = File(outputfilename, "a");
  string cmd = format("~/Github/bwa/bwa mem -v 2 -t 12 -a %s %s", referencepath, fastqpath);
  writefln("Aligning reads from %s to %s using bwa", fastqpath, baseName(referencepath));
  auto ret = executeShell(cmd);
  writefln("Processing %.2f Megabyte of BWA output", to!float(ret.output.length) / (1024 * 1024));
  size_t lines = 0;
  foreach (line; ret.output.splitLines()) {
    if (line[0] == '@') continue;
    lines++;
    auto sline = line.split("\t");
    if (sline.length > 4) {
      if (sline[2] == "*") continue; // Not aligned, just continue with the next line
      string qname = sline[0];
      string rname = sline[2];
      size_t pos = to!size_t(sline[3]);
      size_t mapq = to!size_t(sline[4]);
      if (fq.reads[qname].nAlignments > 0 || mapq > minMapQ) { // If we already had an alignment of the read, the mapQ doesn't matter anymore
        fq.reads[qname].nAlignments++;// Read has been aligned
        fq.reads[qname].alignLine = line;// Read has been aligned
      }
    }
  }
  writefln("Processed %d lines of BWA output", lines);
  size_t mapped = 0;
  auto keys = fq.reads.keys; // Go through the keys and remove the ones that have been uniquely aligned
  foreach (key; keys) {
    if (fq.reads[key].nAlignments == 1) {
      ofp.write(fq.reads[key].alignLine ~ "\n");
      fq.reads.remove(key);
      fq.nUnmapped = fq.nUnmapped - 1;
      mapped++;
    } else { // none or multiple alignments found, reset the nAlignments to 0
      fq.reads[key].nAlignments = 0;
    }
  }
  writefln("Parsed %s lines, mapped %s reads, unmapped reads left %s", lines, mapped, fq.nUnmapped);
}

int main (string[] args) {
  string fastqpath = "/halde/Hi-C/Human/ENCFF319AST.1Mio.fastq.gz";
  string referencepath = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa.gz";
  size_t readLength = 25;

  string outputfilename = "ReadAlignments.txt";
  if(outputfilename.exists) {
    writefln("Deleting previous output file: %s", outputfilename);
    outputfilename.remove();
  }
  FastQ fq = infoFastQ(fastqpath);

  while (readLength <= fq.readLength) {
    if(fq.nUnmapped == 0) break;
    string path = fq.reduceFastQ(readLength);
    fq.mapToGenome(path, referencepath, outputfilename);
    if (path.exists) {
      writefln("Deleting temporary input file: %s", path);
      path.remove();
    }
    readLength += 5;
  }

  return(0);
}

// Iterative fragment mapping 

import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip, zip;
import iopipe.bufpipe : bufd;

import std.io : IOFile = File;
import std.stdio : File, writeln, writefln;
import std.typecons : refCounted;
import std.string : format, split, strip, toUpper, indexOf, splitLines;
import std.conv : to;
import std.zlib : Compress, compress, HeaderFormat, Z_SYNC_FLUSH, Z_FULL_FLUSH;
import std.file : exists, remove;
import std.process : executeShell;
import std.path : baseName;

import containers.hashmap;

//@HISEQ-2500-1:44:C5UH7ANXX:7:1101:15031:68853 2:N:0:

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
  HashMap!(string, MapInfo) reads;
}

//bool has(const HashMap!(string, MapInfo) aa, string k){ return(); }

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
  string outputbuffer;
  foreach (line; fp.assumeText.byLineRange!true) {
    if (nlines % 4 == 0) {
      string shortname = (to!string(line[1 .. $])).split(" ")[0];
      if ((shortname in fq.reads) !is null) { // read in the AA, we need to align it
        notAlignedYet = true;
      } else {
        notAlignedYet = false;
      }
    }
    if (nlines % 4 == 1 || nlines % 4 == 3) { // Read and Quality score need to be reduced
      if (notAlignedYet){
        outputbuffer ~= line[0..readLength] ~ "\n";
      }
    } else {
      if (notAlignedYet){
        outputbuffer ~= line;
      }
    }
    nlines++;
    if(outputbuffer.length > 1024 * 1024){
      ofp.rawWrite(cmp.compress(outputbuffer));
      ofp.rawWrite(cmp.flush(Z_FULL_FLUSH));
      outputbuffer = "";
    }
  }
  ofp.rawWrite(cmp.compress(outputbuffer));
  ofp.rawWrite(cmp.flush());
  return(outputfilename);
}

// Map a fastq file to the reference genome and update the number of unmapped reads (nUnmapped)
void mapToGenome(ref FastQ fq, string fastqpath, string referencepath, string outputfilename = "alignments.txt", size_t minMapQ = 30) {
  auto ofp = File(outputfilename, "a");
  string cmd = format("~/Github/bwa/bwa mem -v 2 -t 12 -T 10 %s %s", referencepath, fastqpath);
  writefln("Aligning reads from %s to %s using bwa", fastqpath, baseName(referencepath));
  auto ret = executeShell(cmd);
  writefln("Processing %.2f Megabyte of BWA output", to!float(ret.output.length) / (1024 * 1024));
  size_t lines = 0;
  size_t unknownreads = 0;
  foreach (line; ret.output.splitLines()) {
    if (line[0] == '@') continue;
    lines++;
    auto sline = line.split("\t");
    if (sline.length > 4) {
      if (sline[2] == "*") continue; // Not aligned, just continue with the next line
      string qname = sline[0];
      if ((qname in fq.reads) is null) {
        unknownreads++;
        continue;
      }
      size_t mapq = to!size_t(sline[4]);
      if (fq.reads[qname].nAlignments > 0 || mapq > minMapQ) { // If we already had an alignment of the read, the mapQ doesn't matter anymore
        fq.reads[qname] = MapInfo(fq.reads[qname].nAlignments + 1, line.idup);
      }
    }
  }
  writefln("Processed %d lines of BWA output (%s unknown reads)", lines, unknownreads);
  //Step 2: figure out which keys need to be removed since they have 1 unique alignment
  size_t mapped = 0;
  unknownreads = 0;
  string[] toremove;
  foreach (key; fq.reads.byKey) {
    if ((key in fq.reads) is null) {
      unknownreads++;
      continue;
    }
    if (fq.reads[key].nAlignments == 1) {
      ofp.write(fq.reads[key].alignLine ~ "\n");
      toremove ~= key;
      mapped++;
    } else { // none or multiple alignments found, reset the nAlignments to 0
      fq.reads[key] = MapInfo(0);
    }
  }

  //Step 3: Remove keys from the AA
  unknownreads = 0;
  foreach(key; toremove){
    if ((key in fq.reads) is null) {
      unknownreads++;
      continue;
    }
    fq.reads.remove(key);
    fq.nUnmapped = fq.nUnmapped - 1;
  }

  writefln("Parsed %s lines, mapped %s reads, unmapped reads left %s", lines, mapped, fq.nUnmapped);
}

// dub -- /halde/Hi-C/Human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /halde/Hi-C/Human/ENCFF319AST.fastq.gz ENCFF319AST.alignment
// dub -- /halde/Hi-C/Human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /halde/Hi-C/Human/ENCFF478EAB.fastq.gz ENCFF478EAB.alignment
int main (string[] args) {
  if(args.length < 4){
    writeln("Please provide the fastq input file and fasta reference");
    return(-1);
  }
  string referencepath = args[1]; //"/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa.gz";
  string fastqpath = args[2]; //"/halde/Hi-C/Human/ENCFF319AST.1Mio.fastq.gz";
  string outputfilename = args[3]; //"ReadAlignments.txt";
  string tmpfmt = baseName(args[3], ".alignment") ~ "%s.fq.gz";
  writefln("tmpfmt: %s", tmpfmt);
  size_t readLength = 25;
  if(outputfilename.exists) {
    writefln("Deleting previous output file: %s", outputfilename);
    outputfilename.remove();
  }
  FastQ fq = infoFastQ(fastqpath);

  while (readLength <= fq.readLength) {
    if(fq.nUnmapped == 0) break;
    string path = fq.reduceFastQ(tmpfmt, readLength);
    fq.mapToGenome(path, referencepath, outputfilename);
    if (path.exists) {
      writefln("Deleting temporary input file: %s", path);
      path.remove();
    }
    readLength += 5;
  }

  return(0);
}

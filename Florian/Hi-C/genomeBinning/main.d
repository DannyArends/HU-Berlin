// Assign paired-end reads that have been mapped to a digested genome to genomic bins
// (c) Danny Arends (HU-Berlin) Sept - 2019
// Run on hulgf34:
// dub -- /halde/Hi-C/Human/genome_bins.txt merged.alignments.txt binned.alignments.txt

import std.stdio : File, writeln, writefln;
import std.conv : to;
import std.string : format, split, strip, toUpper, indexOf, splitLines;

struct Bin {
  size_t start;
  size_t end;
}

struct Chromosome {
  size_t sbin;
  Bin[] bins;
}

struct Genome {
  Chromosome[string] chrs;
  size_t nbins;
}

bool has(T)(const T[string] aa, string k){ return((k in aa) !is null); }

// Read the genome bins into memory
Genome readBins(string path) {
  Genome genome;
  auto fp = File(path, "r");
  size_t nbin = 0;
  //writefln("Starting to read bins from %s", path);
  foreach (line; fp.byLine) {
    auto elements = strip(line).split("\t");
    string chr = to!string(elements[0]);
    if (!has(genome.chrs, chr)) { 
      //writefln("Creating new chromosome %s at bin: %d", chr, nbin);
      genome.chrs[chr] = Chromosome(nbin); 
    }
    genome.chrs[chr].bins ~= Bin(to!size_t(elements[1]), to!size_t(elements[2]));
    nbin = nbin + 1;
  }
  genome.nbins = nbin;
  return(genome);
}

// Check if a digestion fragment is within a bin
int inBin(string chr, size_t pos, Genome genome){
  if ((chr in genome.chrs) is null) return(-1);
  size_t index = genome.chrs[chr].sbin;
  foreach (bin; genome.chrs[chr].bins) {
    if (pos >= bin.start && pos <= bin.end) return(to!int(index));
    index++;
  }
  return(-1);
}

// Bin the aligned reads per digestion fragment into the genomic bin size
void binAlignments(string path, Genome genome, File fpo){
  /*size_t[][] alignments;
  alignments.length = genome.nbins;
  for(size_t i = 0; i < genome.nbins; i++) {
    alignments[i].length = genome.nbins;
  }*/

  auto fp = File(path, "r");
  size_t n = 0;
  size_t a = 0;
  foreach (line; fp.byLine) {
    auto elements = line.split("\t");
    n = n + 1;
    if(to!string(elements[2]) == "NA") continue;
    if(to!string(elements[3]) == "NA") continue;
    if(to!string(elements[5]) == "NA") continue;
    if(to!string(elements[6]) == "NA") continue;
    int i1 = inBin(to!string(elements[1]), (to!size_t(elements[2]) + to!size_t(elements[3])) / 2, genome);
    int i2 = inBin(to!string(elements[4]), (to!size_t(elements[5]) + to!size_t(elements[6])) / 2, genome);
    if (i1 > -1 && i2 > -1 ) {
      //alignments[i1][i2]++;
      fpo.writeln(format("%s\t%d\t%d", to!string(elements[0]), i1, i2));
      a = a + 1;
    }
    if(n % 1000000 == 0) writefln("Done %s of %s", a, n);
  }
  //return(alignments);
}

int main (string[] args) {
  if (args.length < 3) {
    writeln("Please provide:\nGenomic bins [arg1]\nmerged alignments [arg2]\noutput file name [arg3]");
    return(-1);
  }
  string binPath = args[1];  //e.g. "/halde/Hi-C/Human/genome_bins.txt";
  string mergedAlignments = args[2];  //e.g. "merged.alignments.txt";
  string outpath = args[3];  //e.g. "binned.alignments.txt";
  auto fp = File(outpath, "w");

  writefln("Loading bins");
  auto genome = binPath.readBins();
  writefln("Matching alignments");
  binAlignments(mergedAlignments, genome, fp);
  /*writefln("Writing bins");
  for(size_t i = 0; i < genome.nbins; i++) {
    for(size_t j = 0; j < genome.nbins; j++) {
      if(j > 0) fp.write("\t");
      fp.write(alignments[i][j]);
    }
    fp.write("\n");
  } */
  return(0);
}

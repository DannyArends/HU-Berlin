// Matching read (fragments) to an in-silico digested genome
// (c) Danny Arends (HU-Berlin) Sept - 2019

import std.stdio : File, writeln, writefln;
import std.conv : to;
import std.string : format, split, strip, toUpper, indexOf, splitLines;

import containers.hashmap : HashMap;

struct GenomeFragment{
  size_t start;
  size_t end;
}

struct Chromosome {
  GenomeFragment[] fragments;
}

struct Alignment {
  string chromosome;
  size_t position;
  size_t length;
}

bool has(T)(const T[string] aa, string k){ return((k in aa)!is null); }

Chromosome[string] readDigestion(string digestionPath) {
  Chromosome[string] res;
  auto fp = File(digestionPath, "r");
  foreach (line; fp.byLine) {
    auto elements = line.split("\t");
    if (elements.length > 3) {
      string chr = to!string(elements[0]);
      if (!has(res, chr)) {
        res[chr] = Chromosome();
      }
      res[chr].fragments ~= GenomeFragment(to!size_t(elements[1]), to!size_t(elements[2]));
    }
  }
  return(res);
}

HashMap!(string, Alignment) readAlignment(string path) {
  HashMap!(string, Alignment) res;
  auto fp = File(path, "r");
  size_t l = 0;
  foreach (line; fp.byLine) {
    auto elements = line.split("\t");
    res[to!string(elements[0]).idup] = Alignment(to!string(elements[2]).idup, to!size_t(elements[3]), to!string(elements[9]).length);
    l++;
    if(l % 1000000 == 0) writefln("Done %s of %s", l, path);
  }
  return(res);
}

int findFragmentIndex(const Chromosome[string] digestion, const Alignment a) {
  size_t index = 0;
  if ((a.chromosome in digestion) is null) return(-1);
  foreach (fragment; digestion[a.chromosome].fragments) {
    if (a.position >= fragment.start && a.position <= fragment.end) return(to!int(index));
    index++;
  }
  return(-1);
  assert(0, format("Alignment not in bound of chromosome %s", a));
}

// dub -- /halde/Hi-C/Human/HindIII_Digested.bed ENCFF319AST.alignment ENCFF478EAB.alignment merged.alignments.txt

int main (string[] args) {
  if (args.length < 5) {
    writeln("Please provide:\nDigested genome file [arg1]\nAlignment file 1 [arg2]\nAlignment file 2 [arg3]");
    return(-1);
  }
  string digestionPath = args[1];  //e.g. "/halde/Hi-C/Human/HindIII_Digested.bed";
  string path1 = args[2];  //e.g. "ENCFF319AST.alignment";
  string path2 = args[3];  //e.g. "ENCFF478EAB.alignment";
  
  string outfile = args[4];  //e.g. "ENCFF478EAB.alignment";

  auto digestion = readDigestion(digestionPath);
  auto a1 = readAlignment(path1);
  auto a2 = readAlignment(path2);
  foreach (key; digestion.keys) {
    writefln("Chromosome %s has %s digestion fragments", key, digestion[key].fragments.length);
  }
  writefln("Alignment file1: %s aligned reads", a1.keys.length);
  writefln("Alignment file2: %s aligned reads", a2.keys.length);
  string[] paired;
  size_t l = 0;
  size_t p = 0;
  foreach (readname; a1.keys) {
    if (readname in a2) {
      paired ~= readname;
      p++;
    }
    l++;
    if(l % 1000000 == 0) writefln("Matched %s/%s", p, l);
  }
  writefln("Paired a1 in a2: %s", paired.length);
  auto ofp = File(outfile, "w");
  foreach(readname; paired){
    auto chr1 = a1[readname].chromosome;
    auto chr2 = a2[readname].chromosome;
    auto i1 = digestion.findFragmentIndex(a1[readname]);
    string p1 = format("%s\tNA\tNA", chr1);
    if(i1 > 0) p1 = format("%s\t%s\t%s", chr1, digestion[chr1].fragments[i1].start, digestion[chr1].fragments[i1].end);
    auto i2 = digestion.findFragmentIndex(a2[readname]);
    string p2 = format("%s\tNA\tNA", chr2);
    if(i2 > 0) p2 = format("%s\t%s\t%s", chr2, digestion[chr2].fragments[i2].start, digestion[chr2].fragments[i2].end);
    //writefln("%s in chromosomes %s and %s", readname, p1, p2);
    ofp.write(format("%s\t%s\t%s\n", readname, p1, p2));
  }
  return(0);
}

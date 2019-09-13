import std.stdio : File, writeln, writefln;
import std.conv : to;
import std.string : format, split, strip, toUpper, indexOf, splitLines;

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

Alignment[string] readAlignment(string path) {
  Alignment[string] res;
  auto fp = File(path, "r");
  foreach (line; fp.byLine) {
    auto elements = line.split("\t");
    res[to!string(elements[0])] = Alignment(to!string(elements[2]), to!size_t(elements[3]), to!string(elements[3]).length);
  }
  return(res);
}

int findFragmentIndex(const Chromosome[string] digestion, const Alignment a) {
  size_t index = 0;
  foreach(fragment; digestion[a.chromosome].fragments){
    if(a.position >= fragment.start && a.position <= fragment.end) return(to!int(index));
    index++;
  }
  return(-1);
  assert(0, format("Alignment not in bound of chromosome %s", a));
}

// dub -- /halde/Hi-C/Human/HindIII_Digested.bed ENCFF319AST.alignment ENCFF478EAB.alignment

int main (string[] args) {
  if (args.length < 4) {
    writeln("Please provide:\nDigested genome file [arg1]\nAlignment file 1 [arg2]\nAlignment file 2 [arg3]");
    return(-1);
  }
  string digestionPath = args[1];  //e.g. "/halde/Hi-C/Human/HindIII_Digested.bed";
  string path1 = args[2];  //e.g. "ENCFF319AST.alignment";
  string path2 = args[3];  //e.g. "ENCFF478EAB.alignment";
  
  auto digestion = readDigestion(digestionPath);
  auto a1 = readAlignment(path1);
  auto a2 = readAlignment(path2);
  foreach (key; digestion.keys) {
    writefln("Chromosome %s has %s digestion fragments", key, digestion[key].fragments.length);
  }
  writefln("Alignment file1: %s aligned reads", a1.keys.length);
  writefln("Alignment file2: %s aligned reads", a2.keys.length);
  string[] paired;
  foreach(readname; a1.keys){
    if(a2.has(readname)){
      paired ~= readname;
    }
  }
  writefln("Paired a1 in a2: %s", paired.length);
  foreach(readname; paired){
    auto chr1 = a1[readname].chromosome;
    auto chr2 = a2[readname].chromosome;
    auto i1 = digestion.findFragmentIndex(a1[readname]);
    string p1 = format("%s:NA-NA", chr1);
    if(i1 > 0) p1 = format("%s:%s-%s", chr1, digestion[chr1].fragments[i1].start, digestion[chr1].fragments[i1].end);
    auto i2 = digestion.findFragmentIndex(a2[readname]);
    string p2 = format("%s:NA-NA", chr2);
    if(i2 > 0) p2 = format("%s:%s-%s", chr2, digestion[chr2].fragments[i2].start, digestion[chr2].fragments[i2].end);
    writefln("%s in chromosomes %s and %s", readname, p1, p2);
  }
  return(0);
}

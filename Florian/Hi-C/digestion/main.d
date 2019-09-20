// In-silico digestion of a gzip fasta file using a restriction enzyme
// (c) Danny Arends (HU-Berlin) Sept - 2019

import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip;
import iopipe.bufpipe : bufd;

import std.io : File;
import std.stdio : writeln, writefln;
import std.typecons : refCounted;
import std.string : strip, toUpper, indexOf;
import std.conv : to;

struct Enzyme {
  string name;
  string site;
}

enum Enzymes : Enzyme {
  AatII = Enzyme("AatII", "GACGTC"),
  AbsI = Enzyme("AbsI", "CCTCGAGG"),
  Acc65I = Enzyme("Acc65I", "GGTACC"),
  BamHI = Enzyme("BamHI", "GGATCC"),
  EcoRI = Enzyme("EcoRI", "GAATTC"),
  HindIII = Enzyme("HindIII", "AAGCTT"),
  XmnI = Enzyme("XmnI", "GAANNNNTTC")
}

@nogc pure bool compatible (const char a, const char b) nothrow {
  if ("ATGC".indexOf(b) == -1) return false;
  switch(a) {
    case 'A': return b==a;
    case 'T': return b==a;
    case 'G': return b==a;
    case 'C': return b==a;

    case 'W': return "AT".indexOf(b)!=-1;
    case 'S': return "GC".indexOf(b)!=-1;
    case 'R': return "AG".indexOf(b)!=-1;
    case 'Y': return "CT".indexOf(b)!=-1;
    case 'M': return "AC".indexOf(b)!=-1;
    case 'K': return "GT".indexOf(b)!=-1;


    case 'B': return "CGT".indexOf(b)!=-1;
    case 'D': return "AGT".indexOf(b)!=-1;
    case 'H': return "ACT".indexOf(b)!=-1;
    case 'V': return "ACG".indexOf(b)!=-1;

    case 'N': return "ACGT".indexOf(b)!=-1;
    default: return false;
  }
}

void digest (string referencepath, Enzyme enzyme) {
  string chromosome; // Current chromosome we're analyzing
  size_t cPos = 0; // Current position on the chromosome
  size_t pPos = 0; // Previous position of a cut-site
  size_t bLength = 0; // Length of the current buffer
  size_t siteSize = enzyme.site.length; // Length of the cut site
  char[] buffer = new char[](siteSize); // 'Ring' type buffer holding the current genome slice

  auto fp = File(referencepath).refCounted.bufd.unzip(CompressionFormat.gzip);

  foreach (line; fp.assumeText.byLineRange) {
    if (line[0] == '>') { // Fasta header line
      if (cPos > 0) { // Output the fragment that occurs from the last cut-site to the end of chromosome
        writefln("%s\t%s\t%s\t%s\t%s", chromosome, pPos, cPos, enzyme.name, cPos-pPos);
      }
      if (line[1..4] == "CHR") { // Only do the normal chromosomes not the patches
        break;
      } else {
        chromosome = to!string(strip(line[1 .. 3])); // Get the new chromosome name
        cPos = 0; // Chromosome position is 0
        pPos = 0; // Previous cut site is 0
        bLength = 0; // Reset the buffer length
      }
    } else { // Line should contain sequence
      foreach (char c; line) {
        cPos++;
        buffer[bLength] = c;
        bLength++;
        if (bLength == siteSize) { //Buffer is full, compare the sequence
          size_t i = 0;
          for (i = 0; i < buffer.length; i++) {
            if(!compatible(to!char(enzyme.site[i].toUpper()), to!char(buffer[i].toUpper()))) break;
          }
          if (i == buffer.length) { // The whole buffer matched the cut site
            writefln("%s\t%s\t%s\t%s\t%s", chromosome, pPos, cPos - siteSize, enzyme.name, ((cPos - siteSize)-pPos));
            pPos = cPos-siteSize;
          }
          buffer = buffer[1 .. siteSize] ~ "X";
          bLength--;
        }
      }
    }
  }
}

int main (string[] args) {
  digest("/halde/Hi-C/Human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz", Enzymes.HindIII);
  return(0);
}

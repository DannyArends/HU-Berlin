// Analysis of BFMI origin using MGP data
// (c) Danny Arends (HU-Berlin) Sept - 2019

import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip;
import iopipe.bufpipe : bufd;

import std.io : File;
import std.stdio : writeln, writefln;
import std.typecons : refCounted;
import std.string : strip, toUpper, indexOf;
import std.conv : to;

// Call SNPs on the BFMI: nohup ./callSNPs.sh &
int main (string[] args) {
  size_t n = 0;
  auto fp = File("/halde/BFMIorigin/mgp.v5.merged.snps_all.dbSNP142.vcf.gz").refCounted.bufd.unzip(CompressionFormat.gzip);
  foreach (line; fp.assumeText.byLineRange) {
    writefln("%s", line);
    if(n > 300) break;
    n++;
  }
  return(0);
}

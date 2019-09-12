// Iterative fragment mapping 

import iopipe.textpipe : assumeText, byLineRange;
import iopipe.zip : CompressionFormat, unzip;
import iopipe.bufpipe : bufd;

import std.io : File;
import std.stdio : writeln, writefln;
import std.typecons : refCounted;
import std.string : strip, toUpper, indexOf;
import std.conv : to;
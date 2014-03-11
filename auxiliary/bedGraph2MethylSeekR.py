#!/usr/bin/python
import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Convert a series of bedGraph files into input files appropriate for MethylSeekR. Note that the single-base bedGraphs must be first merged with bison_merge_CpGs!')
parser.add_argument('files', metavar='file.bedGraph', nargs='*', help="Input bedGraph file(s).")
args = parser.parse_args()

if((args.files == None) or len(args.files) == 0) :
    parser.print_help()
    sys.exit()

for fname in args.files :
    f = csv.reader(open(fname, "r"), dialect="excel-tab")
    oname=fname.sub(".bedGraph", ".MethylSeekR")
    of = open(oname, "w")
    sys.stdout.write("Processing %s, output will be written to %s\n" % (f, oname))
    sys.stdout.flush()

    for line in f :
        if(line.startswith('track')) :
            continue
        if(int(line[2]) - int(line[1]) != 2) :
            sys.stdout.write("MethylSeekR expects per-CpG metrics and you have per-cytosine metrics! Please run bison_merge_CpGs first and then use the resulting merged bedGraph file(s) as input.\n")
            sys.exit()

        of.write("%s\t%i\t%i\t%i\n" % (line[0], line[1], int(line[1])+1, int(line[4]) + int(line[5]), int(line[4])))
    of.close()

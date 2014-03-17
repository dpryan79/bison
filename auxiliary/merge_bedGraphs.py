#!/usr/bin/python
import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Merge a number of bedGraph files from the bison methylation extractor')
parser.add_argument('outfile', metavar='outfile', help="Output bedGraph files")
parser.add_argument('files', metavar='files', nargs='*', help="Input bedGraph files. There must be at least 2.")
args = parser.parse_args()

if((args.outfile == None) or (args.files == None) or (len(args.files) < 2)) :
    parser.print_help()
    sys.exit()

files = []
for f in args.files :
    if(f != args.outfile) :
        files.append(csv.reader(open(f, "r"), dialect="excel-tab"))
of = open(args.outfile, "w")
of.write("track type=bedGraph\n")

lines = []
for f in files :
    line = f.next()
    if(line[0].startswith('track')) :
        line = f.next()
    lines.append([line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])])

n_finished = 0
n_total = len(files)
while(n_finished < n_total) :
    i = 0
    lowest = 0
    #Determine the appropriate starting point
    while(i<n_total) :
        if(lines[i][0] != None) :
            if(lines[i][0] < lines[lowest][0]) :
                lowest = i
            elif(lines[i][0] == lines[lowest][0] and lines[i][1] < lines[lowest][1]) :
                lowest = i
        elif(lowest == i) :
                lowest += 1
        i += 1

    current = lines[lowest]
    if(lines[lowest][0] == None) :
        print("Oh shit, this shouldn't happen!")
        print(lowest, lines)
        break

    i = 0
    while(i < n_total) :
        if(i != lowest and lines[i][0] != None) :
            if(lines[i][0] == current[0] and lines[i][1] == current[1]) :
                current[3] += lines[i][3]
                current[4] += lines[i][4]
                try :
                    line = files[i].next()
                    lines[i] = [line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])]
                except :
                    lines[i][0] = None
                    n_finished += 1
        i += 1
    frac = round(1000*float(current[3])/float(current[3]+current[4]))
    of.write("%s\t%i\t%i\t%i\t%i\t%i\n" % (current[0], current[1], current[2], frac, current[3], current[4]))

    #Don't forget to increment the "lowest" one too
    try :
        line = files[lowest].next()
        lines[lowest] = [line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])]
    except :
        lines[lowest][0] = None
        n_finished += 1

of.close()

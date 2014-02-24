#!/usr/bin/python
import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Convert a series of bedGraph files into input files appropriate for BSseq.')
parser.add_argument('-chr', metavar='chromosome', help="Only output this chromosome (e.g. chr17) instead of all of them.")
parser.add_argument('prefix', metavar='prefix', help="Output prefix")
parser.add_argument('files', metavar='files', nargs='*', help="Input bedGraph files. There must be at least 2.")
args = parser.parse_args()

if((args.prefix == None) or (args.files == None) or (len(args.files) < 2)) :
    parser.print_help()
    sys.exit()

files = []
for f in args.files :
    files.append(csv.reader(open(f, "r"), dialect="excel-tab"))
ofM = open("%s.M" % (args.prefix), "w")
ofCov = open("%s.Cov" % (args.prefix), "w")
ofbed = open("%s.bed" % (args.prefix), "w")

lines = []
for f in files :
    line = f.next()
    lines.append([line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])])

#Add a header
first = 1
for f in args.files :
    if(first == 1) :
        ofM.write("%s" % f)
        ofCov.write("%s" % f)
        first = 0
    else :
        ofM.write("\t%s" % f)
        ofCov.write("\t%s" % f)
ofM.write("\n")
ofCov.write("\n")

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

    output = 0
    if(args.chr != None) :
        if(lines[lowest][0] == args.chr) :
            ofbed.write("%s\t%i\t%i\n" % (lines[lowest][0], lines[lowest][1], lines[lowest][2])) #Now 1-based
            output = 1
    else :
        ofbed.write("%s\t%i\t%i\n" % (lines[lowest][0], lines[lowest][1], lines[lowest][2])) #Now 1-based
        output = 1

    if(output == 1) :
        i = 0
        while(i < n_total) :
            if(i != 0 and output == 1) :
                ofM.write("\t")
                ofCov.write("\t")
    
            if(lines[i][0] != None) :
                if(lines[i][0] == current[0] and lines[i][1] == current[1]) :
                    ofM.write("%i" % (lines[i][3]))
                    ofCov.write("%i" % (lines[i][3] + lines[i][4]))
                    try :
                        line = files[i].next()
                        lines[i] = [line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])]
                    except :
                        lines[i][0] = None
                        n_finished += 1
                else :
                    ofM.write("0")
                    ofCov.write("0")
            else :
                ofM.write("0")
                ofCov.write("0")

            i += 1
        ofM.write("\n")
        ofCov.write("\n")
    else :
        #We're on the wrong chromosome
        i = 0
        while(i < n_total) :
            if(lines[i][0] != None) :
                if(lines[i][0] == current[0] and lines[i][1] == current[1]) :
                    try :
                        line = files[i].next()
                        lines[i] = [line[0],int(line[1]), int(line[2]), int(line[4]), int(line[5])]
                    except :
                        lines[i][0] = None
                        n_finished += 1
            i += 1

ofM.close()
ofCov.close()
ofbed.close()

##########################################################################################
#                                                                                        #
# Copyright 2018 Megaprobe-Lab                                                           #
#                                                                                        #
# This is software created by the megaprobe lab under the GPL3 license.                  #
#                                                                                        #
# This program parses any ole GFA file you may have with read information for two        #
# organisms. To run the program, utilize the command: "Python2.7 Second_Program.py GFA". #
# The GFA file utilized must be located on the same directory as the parser or full path #
# must be specified.																     #
#                                                                                        #
##########################################################################################

# Load Modules
import os
import sys
import fileinput
import argparse
import gfapy


# Make output directory
if not os.path.exists("Parser_Output"):
    os.mkdir("Parser_Output")


############################################
#                                          #
# Readin GFA format input                  #
#                                          #
############################################


# def coverageSegmentF(C,kmerA,kmerB,Db,k):
# 	global g
# 	kA,kB = kmerA,kmerB

# 	kA = kA.split("A:")
# 	kA = kA[1]
# 	kA = kA.split(',B:')
# 	kA[1] = kA[1].split(")")[0]

#  	kB = kB.split("A:")
# 	kB = kB[1]
# 	kB = kB.split(',B:')
# 	kB[1] = kB[1].split(")")[0]

# 	coverageA = abs(int(kA[0]) - int(kA[1])) 
# 	coverageB = abs(int(kB[0]) - int(kB[1]))
# 	DifEx = (coverageA + coverageB)/2

# 	if DifEx < C:
# 		return False
# 	else:
# 		a = ("L\t%s\t+\t%s\t%s\t%dM\tKC:i:%d"%(kmerA,kmerB,dB,k-1,int(DifEx)))
# 		return a


def positive(C):
	C = int(C)
	if C <= 0:
		sys.exit("Coverage must be positive integer")
	return C



def main():
    global C,g,args
    #maybe clean all the unnatached segments
    output = gfapy.Gfa()
    output.add_line(g.header)
    all_segments = [output.add_line(str(line)) for line in g.segments]
    links = [str(line) for line in g.edges]

    for x in range(len(links)-1):
        link = links[x]
        kmerA = link.split("\t")
        k = kmerA[5]
        dB = kmerA[4]
        kmerB = kmerA[3]
        kmerA = kmerA[1]
        kA,kB = kmerA,kmerB
        kA = kA.split("A:")
        kA = kA[1]
        kA = kA.split(',B:')
        kA[1] = kA[1].split(")")[0]
        kB = kB.split("A:")
        kB = kB[1]
        kB = kB.split(',B:')
        kB[1] = kB[1].split(")")[0]
        coverageA = abs(int(kA[0]) - int(kA[1])) 
        coverageB = abs(int(kB[0]) - int(kB[1]))
        DifEx = (coverageA + coverageB)/2
        if DifEx < C:
            pass
        else:
            output.add_line("L\t%s\t+\t%s\t%s\t%s\tKC:i:%d"%(kmerA,kmerB,dB,k,int(DifEx)))

    filename = os.path.join("Parser_Output",args.output)
    output.to_file(filename)






parser = argparse.ArgumentParser(description="Creates a GFA file with one or two organisms given a kmer size. \
 If coverage is give eliminates below that coverage ")
parser.add_argument("-f", required=True, help="the output file of dbg.py")
parser.add_argument("-output",required=False,default='output.gfa',help="Output GFA file name")
parser.add_argument("-C", required=True,type=positive,help="Eliminates links with Diferential Expresion below coverage")
args = parser.parse_args()
g = gfapy.Gfa.from_file(args.f) 
C = args.C
# To add more organisms add this parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
# change the name and do another call to build and do multiple merge_dicts calls


if __name__ == "__main__":
    main()
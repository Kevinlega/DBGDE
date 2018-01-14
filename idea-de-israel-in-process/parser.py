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
    global args
    segments = []
    newsegments = []
    links = []
    comments = {}
    k = args.k
    linkdict = {}
    C = args.C
    with open(args.f,'r') as file:
        for line in file.readlines():
            line = line.replace("\n","")
            line = line.split("\t")

            if line[0] == 'S':
                segments.append(tuple(line))

            elif line[0] == 'L':
                links.append(tuple(line))
            elif line[0] == '#':
                comments[line[1]] = line[2]

    for x,j in enumerate(links):
        linkdict[x] = tuple([j[1],j[3]])

    for j in range(len(segments)):
        kmer = segments[j][2]
        while x > len(kmer)-k :
            km1 = kmer[x:x+k]
            km2 = kmer[x+1:(x+1)+k]
            c1 = comments[km1]
            c2 = comments[km2]


            c1 = c1.replace("(","")
            c1 = c1.replace(",", " ")
            c1 = c1.replace(")","")
            c1 = c1.split()

            c2 = c2.replace("(","")
            c2 = c2.replace(",", " ")
            c2 = c2.replace(")","")
            c2 = c2.split()

            c1 = abs(int(c1[0])- int(c1[1]))
            c2 = abs(int(c2[0])- int(c2[1]))
            coverage = int((c1 + c2)/2)

            if coverage < C:
                newsegments.append(kmer[0:x+k])
                print("if:%s"%(kmer[0:x+k]))
                kmer = kmer[x+1:]

                x = 0
            else:
                x +=1
                if x > len(kmer)-k:
                    print("else:%s"%(kmer))
                    newsegments.append(kmer)








parser = argparse.ArgumentParser(description="Eliminates below the given coverage.")
parser.add_argument("-k", type=int, required=True, help="the kmer size")
parser.add_argument("-f", required=True, help="the output file of dbg.py")
parser.add_argument("-output",required=False,help="Output GFA file name")
parser.add_argument("-C", required=True,default=0,type=positive,help="Eliminates links with Diferential Expresion below coverage")

args = parser.parse_args()



# To add more organisms add this parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
# change the name and do another call to build and do multiple merge_dicts calls


if __name__ == "__main__":
    main()
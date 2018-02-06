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

def positive(C):
	C = int(C)
	if C <= 0:
		sys.exit("Coverage must be positive integer")
	return C



def main():
    global C,args
    #maybe clean all the unnatached segments
    segments_link = {}
    dummy = os.path.join("Parser_Output",'dummy.gfa')
    filename = os.path.join("Parser_Output",args.output)
    with open(dummy,'w') as output:
        with open(args.f,'r') as file:
            for line in file:
                line_split = line.split()
                if line_split[0] == 'L':

                    kA = line_split[1].split("A:")
                    kA = kA[1]
                    kA = kA.split(',B:')
                    kA[1] = kA[1].split(")")[0]
                    kB = line_split[3].split("A:")
                    kB = kB[1]
                    kB = kB.split(',B:')
                    kB[1] = kB[1].split(")")[0]
                    coverageA = abs(int(kA[0]) - int(kA[1])) 
                    coverageB = abs(int(kB[0]) - int(kB[1]))
                    DifEx = (coverageA + coverageB)/2
                    if DifEx >= C:
                        segments_link[line_split[1]] = None
                        segments_link[line_split[3]] = None
                        output.write("L\t%s\t%s\t%s\t%s\t%s\tKC:i:%d\n"%(line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],int(DifEx)))

                elif line_split[0] == 'S':
                    output.write(line)
                elif line_split[0] == 'H':
                    output.write(line)

    with open(filename,'w') as output:
        with open(dummy,'r') as file:
            for line in file:
                line_split = line.split()
                if line_split[0] == 'S':
                    if line_split[1] in segments_link:
                        output.write(line)
                elif line_split[0] == 'L':
                    output.write(line)
                elif line_split[0] == 'H':
                    output.write(line)

    del segments_link
    os.remove(dummy)


parser = argparse.ArgumentParser(description="Creates a GFA file with one or two organisms given a kmer size. \
 If coverage is give eliminates below that coverage ")
parser.add_argument("-f", required=True, help="the output file of dbg.py")
parser.add_argument("-output",required=False,default='output.gfa',help="Output GFA file name")
parser.add_argument("-C", required=True,type=positive,help="Eliminates links with Diferential Expresion below coverage")
args = parser.parse_args()
C = args.C
# To add more organisms add this parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
# change the name and do another call to build and do multiple merge_dicts calls


if __name__ == "__main__":
    main()
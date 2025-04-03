import os
import re
import argparse
import sys

# based on the faster header given, we can find the id of that in the readlines list
# take the id of the next fasta header as well then that will return the sequence in between

vars = argparse.ArgumentParser(description='Convert gff file to bed format')

vars.add_argument('--gff',dest='gff', type=str, required=True, help='path to gff file')

vars.add_argument('-o','--output',dest='output', type=str, required=True, help='name of the output file')

vars.add_argument('-d','--directory',dest='dir', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

# if an actual output directory is given
if args.dir:

    outdir = args.dir # set output dir

    # check if dir given exists, if not create one
    if not os.path.isdir(f"{outdir}"):
        os.makedirs(f"{outdir}")

# open fasta file (we dont need to close it later)
with open(args.gff,'r') as GffFile, open(f"{outdir}/{args.output}",'w') as output_file:

    # store the fasta data in a list but do not remove new lines
    GffData = GffFile.readlines()

    # loop through Gff data and parse out the info we need for each line
    for L in GffData:

        Line = L.split()

        Chr = Line[0]

        Start = Line[3]
        Stop = Line[4]

        Gene_ID = Line[8].split(';')[0].replace('ID=','')

        output_file.write(f'{Chr}\t{Start}\t{Stop}\t{Gene_ID}\n')
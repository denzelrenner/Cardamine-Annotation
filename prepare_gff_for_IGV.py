import os
import argparse
import logging
import pandas as pd
import sys

vars = argparse.ArgumentParser(description='Output gff for IGV with colours assigned to each locus type')

vars.add_argument('-o','--output_prefix',dest='output_prefix', type=str, required=True, help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=False, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('--gff',dest='input_gff', type=str, required=True, help='gff file produced from barrnap')

vars.add_argument('--locus_colour',dest='locus_colour', nargs='+', type=str, required=False,default=['5_8S','#ff7f0e','18S','#1f77b4','28S','#2ca02c'], help='enter the name of the rDNA locus, and then the hex code to use. e.g 5_8S #FFFFFF 18S #FFFFFF')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# configure the logger we will be using
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{outdir}/prepare_gff_for_IGV.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
Mainlogger = logging.getLogger(__name__)

# hex to rgb function
def hex_to_rgb(hexcode:str):

    return tuple(int(hexcode.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))

# create dict to store locus and its RGB
Locus_RGB = {Locus:hex_to_rgb(Hex) for Locus, Hex in dict(zip(args.locus_colour[0::2],args.locus_colour[1::2])).items()}
    
# create rgb from the hex code
# open fasta file (we dont need to close it later)
with open(args.input_gff,'r') as GffFile, open(f"{outdir}/{args.output_prefix}.gff3",'w') as output_file, \
    open(f"{outdir}/{args.output_prefix}_5_8S.gff3",'w') as output_file_5_8S, open(f"{outdir}/{args.output_prefix}_18S.gff3",'w') as output_file_18S, \
    open(f"{outdir}/{args.output_prefix}_28S.gff3",'w') as output_file_28S:

    # store the fasta data in a list but do not remove new lines
    GffData = GffFile.readlines()

    # loop through Gff data and parse out the info we need for each line
    for Line in GffData:

        Columns = Line.split('\t')

        Info = Columns[8].rstrip('\n')

        # get rDNA locus in line
        Locus = Info.split(';')[0].lstrip('Name=').rstrip('_rRNA')
        RGB = Locus_RGB[Locus]

        # add to info
        # Info += f";color={','.join([str(Val) for Val in RGB])}"
        Info = f"color={','.join([str(Val) for Val in RGB])}"

        # set new columns value
        Columns[8] = Info

        output_file.write(f'{"\t".join(Columns)}\n')

        if Locus == '18S':
            output_file_18S.write(f'{"\t".join(Columns)}\n')
            
        elif Locus == '5_8S':
            output_file_5_8S.write(f'{"\t".join(Columns)}\n')
            
        elif Locus == '28S':
            output_file_28S.write(f'{"\t".join(Columns)}\n')

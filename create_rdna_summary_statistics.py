import os
import argparse
import logging
import pandas as pd
import sys

# This script will take in the fasta file from barrnap as input and identify unique sequences
# This script should take in the fasta file we will subset, as well as the start and stop of the BOI region
# Input
# fasta file of rDNA hits from Barrnap
# output dir and prefix
# start and end coors for BOI region
# chromsome of interest
# query loci

# Output
# a tsv of ID ,count, and sequence
# a fasta file with header as id and count, and then the sequence

# align the sequences together

# based on the faster header given, we can find the id of that in the readlines list
# take the id of the next fasta header as well then that will return the sequence in between

vars = argparse.ArgumentParser(description='Identify unique sequences from a fasta file')

vars.add_argument('--fasta',dest='assembly_fasta', type=str, required=True, help='assembly fasta file')

vars.add_argument('--barrnap',dest='barrnap_fasta', type=str, required=True, help='fasta file produced from barrnap')

# vars.add_argument('-id','--input_dir',dest='input_drectory', type=str, required=False,default=None, help='name of directory which contains the parsed barrnap data')

vars.add_argument('-o','--output_prefix',dest='output_prefix', type=str, required=True, help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('--NOR_start',dest='NORstart', type=int, required=True, help='start coordinate for NOR')

vars.add_argument('--NOR_end',dest='NORend', type=int, required=True, help='end coordinate for NOR')

vars.add_argument('--NOR_chr',dest='NORchromosome', type=str, required=False,default='RL_1\n', help='chromosome that contains the NOR')

vars.add_argument('--query_loci',dest='query_loci', nargs='+', type=str, required=False,default=['5_8S','18S','28S','5S'], help='give the name of the rDNA loci which you want repeat statistics for')

vars.add_argument('--chromosome_prefix',dest='chromosome_prefix',type=str, required=False,default='RL_',help='prefix for chromosome names')

# vars.add_argument('--chromosome',dest='chromosome', type=str, required=False,default='RL_1', help='specify chromosome of interest')

args = vars.parse_args()


# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# configure the logger we will be using
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{outdir}/find_unique_fasta_sequences_global.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
Mainlogger = logging.getLogger('MainLogger')

# create dict with scaffold name as key
def fasta_to_dict(FastaData:list):

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each header to have a 'group' of sequences
    CurrentSeq = []
    CurrentHeader = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # find if the line begins with a gene
        if Line.startswith('>'): 

            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentHeader.replace('>','')]=CurrentSeq
                CurrentHeader=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentHeader.replace('>','')]=CurrentSeq

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

# open fasta file (we dont need to close it later)
with open(args.assembly_fasta,'r') as AssemblyFastaFile,open(args.barrnap_fasta,'r') as BarrnapFastaFile, \
    pd.ExcelWriter(f'{outdir}/MasterSpreadsheet.xlsx', engine='openpyxl') as ExcelOut:


    # store the fasta data for the assembly in a list but do not remove new lines
    AssemblyFastaData = AssemblyFastaFile.readlines()

    # chromosomes and their sequence will be stored here
    AssemblyFastaRefDict = fasta_to_dict(AssemblyFastaData)

    # store the fasta data for the assembly in a list but do not remove new lines
    BarrnapFastaData = BarrnapFastaFile.readlines()

    # chromosomes and their sequence will be stored here
    BarrnapFastaRefDict = fasta_to_dict(BarrnapFastaData)

    # get the length of all chromosomes in the assembly
    AssemblyStatistics = {}

    for Chr in AssemblyFastaRefDict:

        ChrLen = len(''.join([Seq.rstrip() for Seq in AssemblyFastaRefDict[Chr]]))
        AssemblyStatistics[Chr] = ChrLen
        Mainlogger.info(f'Chromosome {Chr.rstrip()} has length {ChrLen}')

     # extract metrics from each scaffold
    # for Scaf in RefDict:

    #     # scaffold length is the sum of all the lines of the fasta we have taken in 
    #     ScaffoldLength = sum([len(Seq.rstrip('\n')) for Seq in RefDict[Scaf]])

    #     # add to metrics dir
    #     ScaffoldMetricsDf[Scaf.rstrip('\n')] = {'ScaffoldLength':ScaffoldLength}


    # total data... outdf= Locus | total rdna length | average length | % in chromosome | copy number
    MainSummaryTable = {'Locus':[],
                        'Copy Number':[],
                        'Unique Species':[],
                        'Average Length (bp)':[],
                        'Total Length (bp)':[],
                        '% of Genome':[]}
    
    # create output dir
    MainSummaryOutDir = f'{outdir}/MainSummary'

    # check if dir given exists, if not create one
    if not os.path.isdir(MainSummaryOutDir):
        os.makedirs(MainSummaryOutDir)
    
    # go through query loci
    for Locus in args.query_loci:

        # open output files
        with open(f"{MainSummaryOutDir}/{args.output_prefix}_{Locus}.tsv",'w') as TsvFile, \
        open(f"{MainSummaryOutDir}/{args.output_prefix}_{Locus}.fa",'w') as FastaOutput:

            # output lists
            AllFastaSeqs = []

            # loop through ref dict
            for Header in BarrnapFastaRefDict:

                # check if we are dealing with the right query locus
                if Locus in Header:

                    Mainlogger.info(f'Found match: {Header.rstrip()}')

                    # break down the header
                    Start,Stop = Header.split(':')[-1][:-4].split('-')

                    Mainlogger.info(f'Start={Start}, Stop={Stop}\n')

                    # add the sequence to all FastaSeqs
                    AllFastaSeqs.append(BarrnapFastaRefDict[Header][0])

            # get the unique fasta sequences
            UniqueFastaSeqs = list(set(AllFastaSeqs))
            CopyNumber = len(AllFastaSeqs)
            UniqueSpecies = len(UniqueFastaSeqs)
            TotalLength = sum([len(Seq.rstrip()) for Seq in AllFastaSeqs]) if len(AllFastaSeqs) != 0 else 0
            AverageLength = round(TotalLength/CopyNumber,3) if CopyNumber != 0 else 0
            PercentageOfGenome = round(TotalLength/sum(list(AssemblyStatistics.values())),3)

            MainSummaryTable['Locus'].append(Locus)
            MainSummaryTable['Copy Number'].append(CopyNumber)
            MainSummaryTable['Unique Species'].append(UniqueSpecies)
            MainSummaryTable['Average Length (bp)'].append(AverageLength)
            MainSummaryTable['Total Length (bp)'].append(TotalLength)
            MainSummaryTable['% of Genome'].append(PercentageOfGenome)

            # assign a unique ID to each fasta file
            IDFastaSeq = {UniqueFastaSeqs[i]:f'Species_{i}' for i in range(len(UniqueFastaSeqs))}

            # write headers for tsv file
            TsvFile.write(f'ID\tCount\tSequence\n')
        
            # go through each of the unique headers and output their count
            for Seq in UniqueFastaSeqs:

                # get the sequence count
                Frequency = AllFastaSeqs.count(Seq)

                # write id and sequence to tsv
                TsvFile.write(f'{IDFastaSeq[Seq]}\t{Frequency}\t{Seq}')

                # make header for fasta file
                Header = f'>{IDFastaSeq[Seq]} | {Frequency} Sequences Found\n'

                # write to fasta file
                FastaOutput.write(Header)
                FastaOutput.write(Seq)

    # turn chromosome summary table to df
    MainSummaryTableDf = pd.DataFrame(MainSummaryTable)

    MainSummaryTableDf.to_excel(excel_writer=ExcelOut,sheet_name='Genome',index=False)
    MainSummaryTableDf.to_excel(excel_writer=f'{MainSummaryOutDir}/Genome.xlsx',index=False)
    MainSummaryTableDf.to_csv(path_or_buf=f'{MainSummaryOutDir}/Genome.tsv',sep='\t',index=False)






    ### NOR data
    NOROutDir = f'{outdir}/NORSummary'

    # check if dir given exists, if not create one
    if not os.path.isdir(NOROutDir):
        os.makedirs(NOROutDir)

    # total data... outdf= Locus | total rdna length | average length | % in chromosome | copy number
    NORSummaryTable = {'Locus':[],
                        'Copy Number':[],
                        'Unique Species':[],
                        'Average Length (bp)':[],
                        'Total Length (bp)':[],
                        '% of Chromosome':[],
                        '% of Genome':[]}
    
    # go through query loci
    for Locus in args.query_loci:

        # open output files
        with open(f"{NOROutDir}/{args.output_prefix}_{Locus}.tsv",'w') as TsvFile, \
        open(f"{NOROutDir}/{args.output_prefix}_{Locus}.fa",'w') as FastaOutput:

            # output lists
            AllFastaSeqs = []

            # loop through ref dict
            for Header in BarrnapFastaRefDict:

                # if the query locus is not in the BOI, output error?

                # check if we are dealing with the right query locus
                if Locus in Header and args.NORchromosome.rstrip() in Header:

                    Mainlogger.info(f'Found match: {Header.rstrip()}')

                    # break down the header
                    Start,Stop = Header.split(':')[-1][:-4].split('-')

                    Mainlogger.info(f'Start={Start}, Stop={Stop}\n')

                    # # now check if the Start and Stop are within the BOIs
                    if (int(Start) >= args.NORstart and int(Start) <= args.NORend) and (int(Stop) >= args.NORstart and int(Stop) <= args.NORend):

                        # add the sequence to all FastaSeqs
                        AllFastaSeqs.append(BarrnapFastaRefDict[Header][0])

            # get the unique fasta sequences
            UniqueFastaSeqs = list(set(AllFastaSeqs))
            CopyNumber = len(AllFastaSeqs)
            UniqueSpecies = len(UniqueFastaSeqs)
            TotalLength = sum([len(Seq.rstrip()) for Seq in AllFastaSeqs]) if len(AllFastaSeqs) != 0 else 0
            AverageLength = round(TotalLength/CopyNumber,3) if CopyNumber != 0 else 0
            PercentageOfChromosome = round(TotalLength/AssemblyStatistics[args.NORchromosome],3)
            PercentageOfGenome = round(TotalLength/sum(list(AssemblyStatistics.values())),3)

            NORSummaryTable['Locus'].append(Locus)
            NORSummaryTable['Copy Number'].append(CopyNumber)
            NORSummaryTable['Unique Species'].append(UniqueSpecies)
            NORSummaryTable['Average Length (bp)'].append(AverageLength)
            NORSummaryTable['Total Length (bp)'].append(TotalLength)
            NORSummaryTable['% of Chromosome'].append(PercentageOfChromosome)
            NORSummaryTable['% of Genome'].append(PercentageOfGenome)

            # assign a unique ID to each fasta file
            IDFastaSeq = {UniqueFastaSeqs[i]:f'Species_{i}' for i in range(len(UniqueFastaSeqs))}

            # write headers for tsv file
            TsvFile.write(f'ID\tCount\tSequence\n')
        
            # go through each of the unique headers and output their count
            for Seq in UniqueFastaSeqs:

                # get the sequence count
                Frequency = AllFastaSeqs.count(Seq)

                # write id and sequence to tsv
                TsvFile.write(f'{IDFastaSeq[Seq]}\t{Frequency}\t{Seq}')

                # make header for fasta file
                Header = f'>{IDFastaSeq[Seq]} | {Frequency} Sequences Found\n'

                # write to fasta file
                FastaOutput.write(Header)
                FastaOutput.write(Seq)

    # turn chromosome summary table to df
    NORSummaryTableDf = pd.DataFrame(NORSummaryTable)

    NORSummaryTableDf.to_excel(excel_writer=ExcelOut,sheet_name='NOR',index=False)
    NORSummaryTableDf.to_excel(excel_writer=f'{NOROutDir}/NOR.xlsx',index=False)
    NORSummaryTableDf.to_csv(path_or_buf=f'{NOROutDir}/NOR.tsv',sep='\t',index=False)





    ## Per Chromosome Data

    # loop through chromosomes
    for Chr in AssemblyFastaRefDict:

        # create output dir
        ChrDir = f'{outdir}/{Chr.rstrip()}'

        # make output dir
        if not os.path.isdir(ChrDir):
                os.makedirs(ChrDir)

        # create dict for chromosome specific data...  outdf = Locus | total rdna length | average length | % in chromosome | copy number
        ChromosomeSummaryTable = {'Locus':[],
                                'Copy Number':[],
                                'Unique Species':[],
                                'Average Length (bp)':[],
                                'Total Length (bp)':[],
                                '% of Chromosome':[]}
        
        # go through query loci
        for Locus in args.query_loci:

            # create output dir
            QueryLociDir = f'{ChrDir}/{Locus}'

            if not os.path.isdir(QueryLociDir):
                os.makedirs(QueryLociDir)

            # configure the logger we will be using
            logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{QueryLociDir}/find_unique_fasta_sequences_{Locus}.log",mode='w')],datefmt='%H:%M:%S',force=True)

            # create a logger
            logger = logging.getLogger('SubLogger')

            # open output files
            with open(f"{QueryLociDir}/{args.output_prefix}_{Locus}.tsv",'w') as TsvFile, \
            open(f"{QueryLociDir}/{args.output_prefix}_{Locus}.fa",'w') as FastaOutput:

                # output lists
                AllFastaSeqs = []

                # loop through ref dict
                for Header in BarrnapFastaRefDict:

                    # if the query locus is not in the BOI, output error?

                    # check if we are dealing with the right query locus
                    if Locus in Header and Chr.rstrip() in Header:

                        logger.info(f'Found match: {Header.rstrip()}')

                        # break down the header
                        Start,Stop = Header.split(':')[-1][:-4].split('-')

                        logger.info(f'Start={Start}, Stop={Stop}\n')

                        # # now check if the Start and Stop are within the BOIs
                        # if (int(Start) >= args.BOIstart and int(Start) <= args.BOIend) and (int(Stop) >= args.BOIstart and int(Stop) <= args.BOIend):

                        # add the sequence to all FastaSeqs
                        AllFastaSeqs.append(BarrnapFastaRefDict[Header][0])
                        # AllFastaSeqs.append(''.join([Seq.rstrip() for Seq in BarrnapFastaRefDict[Header]]))

                # get the unique fasta sequences
                UniqueFastaSeqs = list(set(AllFastaSeqs))
                CopyNumber = len(AllFastaSeqs)
                UniqueSpecies = len(UniqueFastaSeqs)
                TotalLength = sum([len(Seq.rstrip()) for Seq in AllFastaSeqs]) if len(AllFastaSeqs) != 0 else 0
                AverageLength = round(TotalLength/CopyNumber,3) if CopyNumber != 0 else 0
                PercentageOfChromosome = round(TotalLength/AssemblyStatistics[Chr],3)

                logger.info(f'Number of Total {Locus} sequences in BOI: {len(AllFastaSeqs)}')
                logger.info(f'Number of Unique {Locus} sequences in BOI: {len(UniqueFastaSeqs)}')
                
                # assign a unique ID to each fasta file
                IDFastaSeq = {UniqueFastaSeqs[i]:f'Species_{i}' for i in range(len(UniqueFastaSeqs))}

                # write headers for tsv file
                TsvFile.write(f'ID\tCount\tSequence\n')
            
                # go through each of the unique headers and output their count
                for Seq in UniqueFastaSeqs:

                    # get the sequence count
                    Frequency = AllFastaSeqs.count(Seq)

                    # write id and sequence to tsv
                    TsvFile.write(f'{IDFastaSeq[Seq]}\t{Frequency}\t{Seq}')

                    # make header for fasta file
                    Header = f'>{IDFastaSeq[Seq]} | {Frequency} Sequences Found\n'

                    # write to fasta file
                    FastaOutput.write(Header)
                    FastaOutput.write(Seq)

                logger.info(f'Done')

                # calculate:
                # total length of all rDNA type
                # average length of rDNA type. total length od rDNA type / number of rDNA type
                # percentage of rDNA in genome. total length of rdna/ total genome length
                ChromosomeSummaryTable['Locus'].append(Locus)
                ChromosomeSummaryTable['Copy Number'].append(CopyNumber)
                ChromosomeSummaryTable['Unique Species'].append(UniqueSpecies)
                ChromosomeSummaryTable['Average Length (bp)'].append(AverageLength)
                ChromosomeSummaryTable['Total Length (bp)'].append(TotalLength)
                ChromosomeSummaryTable['% of Chromosome'].append(PercentageOfChromosome)

        # turn chromosome summary table to df
        ChromosomeSummaryTableDf = pd.DataFrame(ChromosomeSummaryTable)

        # output as excel file and tsv
        ChromosomeSummaryTableDf.to_excel(excel_writer=f'{ChrDir}/Chromosome{Chr.rstrip().replace(args.chromosome_prefix,'')}.xlsx',sheet_name=f'Chromosome{Chr.rstrip().replace(args.chromosome_prefix,'')}',index=False)
        ChromosomeSummaryTableDf.to_excel(excel_writer=ExcelOut,sheet_name=f'Chromosome{Chr.rstrip().replace(args.chromosome_prefix,'')}',index=False)
        ChromosomeSummaryTableDf.to_csv(path_or_buf=f'{ChrDir}/Chromosome{Chr.rstrip().replace(args.chromosome_prefix,'')}.tsv',sep='\t',index=False)

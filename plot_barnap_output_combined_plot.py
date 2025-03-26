import os
import re
import argparse
import pandas as pd
import sys


# remotes::install_github("wilkelab/cowplot")
#     install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#     remotes::install_github("clauswilke/colorblindr")

# based on the faster header given, we can find the id of that in the readlines list synteny map, zooming in on BOIs, presence absence gnee profiles of both haplotypes, gene syntney on them, align to herseut

# take the id of the next fasta header as well then that will return the sequence in between

vars = argparse.ArgumentParser(description='Process Barnap Output and Create Plots')

vars.add_argument('--gff',dest='gff', type=str, required=True, help='the gff file produced by barrnap')

vars.add_argument('-fa','--fasta',dest='fasta', type=str, required=True, help='the fasta file containing all the scaffolds')

vars.add_argument('-d','--directory',dest='dir', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('-sw','--sliding_window',dest='sw', type=int, required=False,default=1000, help='size for the sliding window within which metrics are calculated. (default=1000)')

vars.add_argument('-bp','--bp_size',dest='bp_size', type=str, required=False,default='b', help='sliding window value will be interpereted as being in b, kb, or mb. (default=b)')

vars.add_argument('--flip',dest='flip_coors', action='store_true',required=False,default=False, help='set this argument (i.e --flip) to flip the coordinates of the dot plot')

vars.add_argument('--ignore_zero',dest='ignore_zero', action='store_true',required=False,default=False, help='set this argument (i.e --ignore_zero) to stop density values of 0 from being plotted')

vars.add_argument('--query_loci',dest='query_loci', nargs='+', type=str, required=False,default='all_loci', help='give the name of the rDNA loci repeats in order of their appearance in the repeat block (i.e 5s or 5s 28s)')

vars.add_argument('--combined_plot',dest='combined_plot', action='store_true',required=False,default=False, help='set this argument to get a combined plot with eahc scaffold')

args = vars.parse_args()

# print(args.query_loci)
# sys.exit()

# by default directory to write output to is current wd
OutDir = os.getcwd()

if args.dir:
    OutDir = args.dir # set output dir

# define scripts dir and plots dir
ScriptsDirPath = f"{OutDir}/Scripts"
PlotsDirPath = f"{OutDir}/Plots"
OutputDataDirPath = f"{OutDir}/Output_Data"

for Dir in [OutDir,ScriptsDirPath,PlotsDirPath,OutputDataDirPath]:

    # check if dir given exists, if not create one
    if not os.path.isdir(Dir):
        os.makedirs(Dir)

SlidingWindowSize = args.sw

# sort out what the sliding window value is
if args.bp_size != 'b':

    if args.bp_size == 'kb':
        SlidingWindowSize = args.sw * (10**3)
    
    elif args.bp_size == 'mb':
        SlidingWindowSize = args.sw * (10**6)
    
    else:
        sys.exit("Wrong argument provided. -bp only takes b, kb, or mb")


# store the R installation chunk into a variable
# making the R file

# store installation code in its own chunk because it is hard to write
installation_chunk = '''if (!require("ggplot2", quietly = TRUE)) {
    print("The package ggplot2 is not installed. Beginning installation...")
    install.packages("ggplot2",dependencies=TRUE,repos='https://cloud.r-project.org')
    } else {
    print("The package ggplot2 is installed.")
    }

    if (!require("gridExtra", quietly = TRUE)) {
    print("The package gridExtra is not installed. Beginning installation...")
    install.packages("gridExtra",dependencies=TRUE,repos='https://cloud.r-project.org')
    } else {
    print("The package gridExtra is installed.")
    }

    if (!require("ggsci", quietly = TRUE)) {
    print("The package ggsci is not installed. Beginning installation...")
    install.packages("ggsci",repos='https://cloud.r-project.org')
    } else {
    print("The package ggsci is installed.")
    }

    if (!require("colorspace", quietly = TRUE)) {
    print("The package colorspace is not installed. Beginning installation...")
    install.packages("colorspace", repos = 'https://cloud.r-project.org')
    } else {
    print("The package colorspace is installed.")
    }

    if (!require("remotes", quietly = TRUE)) {
    print("The package remotes is not installed. Beginning installation...")
    install.packages("remotes", repos = 'https://cloud.r-project.org')
    } else {
    print("The package remotes is installed.")
    }

    if (!require("ggpubr", quietly = TRUE)) {
    print("The package ggpubr is not installed. Beginning installation...")
    install.packages("ggpubr", repos = 'https://cloud.r-project.org')
    } else {
    print("The package ggpubr is installed.")
    }\n'''

# function which takes in a fasta file and returns a dictionary with fasta header as key and the sequence as the value
def fasta_to_dict(FastaData:list):

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each scaffold to have a 'group' of sequences
    CurrentSeq = []
    CurrentScaffold = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # only run this code if the current gene is 

        # find if the line begins with a gene
        if Line.startswith('>'): 

            
            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentScaffold.replace('>','')]=CurrentSeq
                CurrentScaffold=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentScaffold.replace('>','')]=CurrentSeq

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

# function to grep out the product from the line in 
def grep_product(line:str):

    # grep the pattern
    MatchObj=re.search(r'product.*S',line.split()[0]) # checking for pattern only in first portion of the string
    # MatchObj=re.search(r'product=.*S',line.split()[0])

    # return the specific product found 
    return line[MatchObj.start():MatchObj.end()].replace('product=','').rstrip()

# function to calculate repeat density. returns pandas dataframe with normalised repeat density values
def calculate_repeat_density(input:pd.DataFrame):

    FinalDf = input

    # find the maximum repeat density
    MaxRepeatDensity = max(FinalDf['RepeatDensity'])

    if MaxRepeatDensity == 0:
        MaxRepeatDensity = 1

    # find total number of repeats in the genome
    TotalNumberOfRepeats = sum(FinalDf['RepeatNumber'])

    # prevent 0 dividion error
    if TotalNumberOfRepeats == 0:
        TotalNumberOfRepeats = 1

    # now we can calculate the normalised value
    # NormalisedRepeatDensity = [rd/MaxRepeatDensity if MaxRepeatDensity!=0 else 0 for rd in FinalDf['RepeatDensity']]
    NormalisedRepeatDensity = [rd/TotalNumberOfRepeats for rd in FinalDf['RepeatDensity']]

    # add column to the final df
    FinalDf['NormalisedRepeatDensity'] = NormalisedRepeatDensity

    return FinalDf

# global variable which has Scaffold metrics. format {'Scaffold':{'length(in bp)':200000}}
ScaffoldMetricsDf = {}

# open fasta file (we dont need to close it later)
with open(args.fasta,'r') as FastaFile:

    # store the fasta data in a list but do not remove new lines
    FastaData = FastaFile.readlines()

    # get dict with fasta header as key and sequence as value
    RefDict = fasta_to_dict(FastaData=FastaData)

    # extract metrics from each scaffold
    for Scaf in RefDict:

        # scaffold length is the sum of all the lines of the fasta we have taken in 
        ScaffoldLength = sum([len(Seq.rstrip('\n')) for Seq in RefDict[Scaf]])

        # add to metrics dir
        ScaffoldMetricsDf[Scaf.rstrip('\n')] = {'ScaffoldLength':ScaffoldLength}

# sort order of scaffolds and print out metrics. In future adjust to go with log file
{print(f'Scaffold {k} has a length of {v["ScaffoldLength"]}') for k, v in sorted(ScaffoldMetricsDf.items(), key=lambda Length:Length[1]['ScaffoldLength'],reverse=True)}

print(f"Total Assembly Length:{sum([ScaffoldMetricsDf[chr]['ScaffoldLength'] for chr in ScaffoldMetricsDf])}")

# store data in pandas df
GffOutput = pd.read_csv(filepath_or_buffer=f"{args.gff}",sep='\t',skiprows=1,header=None)

# df stats
# column 0 = scaffold
# column 3 = start
# column 4 = stop
# column 5 = evalue
# column 6 = strand
# column 8 = product

# get the unqiue names for the different scaffolds
ScaffoldNames = GffOutput[0].unique().tolist()

# variable to store different dataframes
CollectiveDfs = {}

# create data frame which we will add data to
# OutputDataFrameAsDict = pd.DataFrame(data={},columns=['Scaffold','GenomicCoordinate','RepeatDensity'])
# OutputDataFrameAsDict = {'Scaffold':[],'GenomicCoordinate':[],'RepeatDensity':[]}

# open plotting file
with open(f"{ScriptsDirPath}/MainCombinedPlot.R",'w') as MainRplotter_file:

    # subset the initial gff to have dfs for each scaffold
    for Scaf in ScaffoldNames:

        # create data frame which we will add data to
        ArryOfOutputDataDicts = []
    
        # create a plotting file for each scaffold separately
        with open(f"{ScriptsDirPath}/Scaf_{Scaf.lstrip('RL_')}_Plot.R",'w') as Rplotter_file:

            TempDf = GffOutput[GffOutput[0] == Scaf]

            # extract only what the rRNA product is in column 8
            ProductColumn = [grep_product(line=ProductLine) for ProductLine in TempDf[8]]

            # drop current column with product info
            FinalDf = TempDf.drop(columns=8)

            # put new column in its place
            FinalDf[8] = ProductColumn

            # first we need to start the sliding window
            SlidingWindowStart = 0
            SlidingWindowEnd = 0

            # print(TempDf)
            # sys.exit()

            # print(FinalDf)

            for GenomicCoor in range(0,ScaffoldMetricsDf[Scaf]["ScaffoldLength"],SlidingWindowSize):

                # get the coors which we will use to subset data from the Df
                SlidingWindowStart = GenomicCoor
                SlidingWindowEnd = SlidingWindowStart + SlidingWindowSize

                # if the sliding window end is larger than the scaffold length we end the loop here and finalise the dat
                if SlidingWindowEnd >= ScaffoldMetricsDf[Scaf]["ScaffoldLength"]:
                    SlidingWindowEnd = ScaffoldMetricsDf[Scaf]["ScaffoldLength"]

                # subset df to get hits in the gff that have a start position within this range
                HitsDf = FinalDf[(FinalDf[3] >= SlidingWindowStart) & (FinalDf[3] < SlidingWindowEnd)]

                # the position for the window 
                RepeatWindowCoor = GenomicCoor

                # repeat density is calculated as the number of hits/repeats in this window divided by the window size
                RepeatDensity = HitsDf.shape[0]/SlidingWindowSize

                # repeat counter variable set here
                RepeatCounter = HitsDf.shape[0]

                # calculate repeat density by finding how many of the query loci are in the window
                if args.query_loci != 'all_loci':

                    # get the expected sequence
                    ExpectedSequence = args.query_loci

                    # set the Sequence Builder
                    SequenceBuilder = []

                    # set a counter for number of repeats
                    RepeatCounter = 0


                    # create dictionary with each loci we are looking for as key
                    # QueryLociCount = {ql:0 for ql in args.query_loci}

                    # grab the column which has the rDNA loci and count how many times we see that repeat 
                    LociInWindow = HitsDf[8]

                    for Loci in LociInWindow:
                        
                        # simple code for quick plots of 5S, 18S etc
                        # if Loci in ExpectedSequence:
                        #     RepeatCounter += 1
                        
                        if Loci in ExpectedSequence:

                            SequenceBuilder.append(Loci)

                            # # check if the position of the loci that we just added is equal to what its position should be
                            if SequenceBuilder.index(Loci) != ExpectedSequence.index(Loci):
                                print(f'Removing...')
                                print(SequenceBuilder)
                                SequenceBuilder = []

                            # check if the sequence builder is expected
                            if SequenceBuilder == ExpectedSequence:
                                RepeatCounter += 1
                                SequenceBuilder = []

                        # if the loci we expect to be at the start is not appended to the start then clear the sequence builder and put it as the first element
                    
                    # now repeat density is the counted loci divided by total window size
                    RepeatDensity = RepeatCounter

                    # now if we are looking for a repeat of more than 1 locus (i.e 5.8,18,28) it is possible that sometimes 
                    
                

                # create row to submit to list of dicts
                RowToSubmit = {'Scaffold':Scaf,
                            'GenomicCoordinate':RepeatWindowCoor,
                            'RepeatDensity':RepeatDensity,
                            'RepeatNumber':RepeatCounter}

                # add to array of all rows
                ArryOfOutputDataDicts.append(RowToSubmit)

                # if SlidingWindowEnd == ScaffoldMetricsDf[Scaf]["ScaffoldLength"]:
                #     print(SlidingWindowStart)
                #     print(SlidingWindowEnd)
                #     print(SlidingWindowSize)
                #     print(HitsDf)
                #     print(RepeatDensity)

                # sys.exit()

            # table containing scaffold name, repeat density and genomic coor of the repeat
            IntermediaryDataFrame = pd.DataFrame(ArryOfOutputDataDicts)
            
            # the final data frame that will be passed as input to R
            OutputDataFrame = calculate_repeat_density(input=IntermediaryDataFrame)


            # add this DF to dict of all DFs
            # CollectiveDfs[Scaf] = FinalDf

            # write the subset df to output file
            # FinalDf.to_csv(path_or_buf=f'{OutputDataDirPath}/Scaf_{Scaf}_Subset.gff',sep='\t',index=False)

            # write table to output
            OutputDataFrame.to_csv(path_or_buf=f'{OutputDataDirPath}/Scaf_{Scaf}_Subset.tsv',sep='\t',index=False)
            
            
            

            

            # write installation code to file
            Rplotter_file.write(installation_chunk)

            # library packages we will be using
            Rplotter_file.write('library(ggplot2)\n')
            Rplotter_file.write('library(gridExtra)\n')
            Rplotter_file.write('library(colorspace)\n')
            Rplotter_file.write('library(colorblindr)\n')
            Rplotter_file.write('library(scales)\n')

            # fill='#69b3a2', color='#e9ecef', alpha=0.8
            # set a working directory
            Rplotter_file.write(f"setwd('{PlotsDirPath}')\n")

            # load in the different files
            Rplotter_file.write(f"df <- read.table(file='{os.path.abspath(OutputDataDirPath)}/Scaf_{Scaf}_Subset.tsv',header=TRUE,sep='\\t')\n")

            # start graphics driver
            # load graphics driver
            Rplotter_file.write(f"pdf('Scaf_{Scaf}_Plot_{args.sw}_{args.bp_size}_Win.pdf', width = 25, height = 25)\n")

            # pos genes plot
            Rplotter_file.write("ggplot(df, aes(x = GenomicCoordinate,y=NormalisedRepeatDensity)) +\n")
            Rplotter_file.write("geom_point(size=4) +\n")
            Rplotter_file.write(f"labs(x = 'Scaffold {Scaf.lstrip('RL_')} position (Mb)',y = 'Normalised Density') +\n")
            
            # flip x and y coors if specified
            if args.flip_coors:
                Rplotter_file.write("coord_flip()+\n")

            # remove grid from background of plot
            # remove grid from background of plot
            Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
            Rplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
            Rplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
            Rplotter_file.write('''theme(axis.title.y = element_text(vjust=-2.5,margin = unit(c(t=0, r=30, b=0, l=0), "mm")),
                axis.title.x = element_text(vjust=5.5,margin = unit(c(t=30, r=0, b=0, l=0), "mm")),axis.text = element_text(size = 35),
                legend.text=element_text(size=40),axis.title=element_text(size = 45),
                legend.title = element_text(size=40),legend.key.size = unit(2, 'cm'),
                axis.ticks.y = element_line(linewidth=2), axis.line.y = element_line(linewidth=2),
                axis.ticks.x = element_line(linewidth=2),axis.line.x = element_line(linewidth=2))+\n''')
            # Rplotter_file.write("scale_fill_manual(labels = c('Pan-genome', 'Positively Correlated Genes'),values=c('darkorange','purple')) +\n")
            # Rplotter_file.write("scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+\n")
            # flip axis
            #Rplotter_file.write("coord_flip()+\n")

            # centre plot tile
            Rplotter_file.write("theme(plot.title=element_text(hjust=0.5))+\n") 

            # make tick marks longer
            Rplotter_file.write("theme(axis.ticks.length=unit(1,'cm'))+\n") 

            # add space to top and right hand side of plot
            Rplotter_file.write("theme(plot.margin=unit(c(t=30,r=30,b=0,l=0),'mm'))+\n")

            # add plot title
            # Rplotter_file.write(f"ggtitle('Placeholder')+\n") 
        
            
            Rplotter_file.write(f"scale_x_continuous(limits=c(0,{ScaffoldMetricsDf[Scaf]['ScaffoldLength']}),labels = label_number(scale= 1e-6))\n") 

            # set limits on x axis
            # Rplotter_file.write(f"scale_x_discrete(limits=c(0, {ScaffoldMetricsDf[Scaf]['ScaffoldLength']}))\n")
                            

            Rplotter_file.write("dev.off()\n")



            # for combined plot code
            # write installation code to file
            MainRplotter_file.write(installation_chunk)

            # library packages we will be using
            MainRplotter_file.write('library(ggplot2)\n')
            MainRplotter_file.write('library(gridExtra)\n')
            MainRplotter_file.write('library(colorspace)\n')
            MainRplotter_file.write('library(colorblindr)\n')
            MainRplotter_file.write('library(scales)\n')
            MainRplotter_file.write('library(ggpubr)\n')
            MainRplotter_file.write('library("gridExtra")\n')
            

            # fill='#69b3a2', color='#e9ecef', alpha=0.8
            # set a working directory
            MainRplotter_file.write(f"setwd('{os.path.abspath(PlotsDirPath)}')\n")

            # load in the different files
            MainRplotter_file.write(f"df_{Scaf.lstrip('RL_')} <- read.table(file='{os.path.abspath(OutputDataDirPath)}/Scaf_{Scaf}_Subset.tsv',header=TRUE,sep='\\t')\n")

            # start graphics driver
            # load graphics driver
            # MainRplotter_file.write(f"pdf('Scaf_{Scaf}_Plot_{args.sw}_{args.bp_size}_Win.pdf', width = 25, height = 25)\n")

            # pos genes plot
            print(f"p_{Scaf.lstrip('RL_')}")
            MainRplotter_file.write(f"p_{Scaf.lstrip('RL_')} <- ggplot(df_{Scaf.lstrip('RL_')}, aes(x = GenomicCoordinate,y=NormalisedRepeatDensity)) +\n")
            MainRplotter_file.write("geom_point(size=4) +\n")
            MainRplotter_file.write(f"labs(x = 'Scaffold {Scaf.lstrip('RL_')} position (Mb)',y = 'Normalised Density') +\n")
            
            # flip x and y coors if specified
            if args.flip_coors:
                MainRplotter_file.write("coord_flip()+\n")

            # remove grid from background of plot
            # remove grid from background of plot
            MainRplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
            MainRplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
            MainRplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
            MainRplotter_file.write('''theme(axis.title.y = element_text(vjust=-2.5,margin = unit(c(t=0, r=30, b=0, l=0), "mm")),
                axis.title.x = element_text(vjust=5.5,margin = unit(c(t=30, r=0, b=0, l=0), "mm")),axis.text = element_text(size = 35),
                legend.text=element_text(size=40),axis.title=element_text(size = 45),
                legend.title = element_text(size=40),legend.key.size = unit(2, 'cm'),
                axis.ticks.y = element_line(linewidth=2), axis.line.y = element_line(linewidth=2),
                axis.ticks.x = element_line(linewidth=2),axis.line.x = element_line(linewidth=2))+\n''')
            # Rplotter_file.write("scale_fill_manual(labels = c('Pan-genome', 'Positively Correlated Genes'),values=c('darkorange','purple')) +\n")
            # Rplotter_file.write("scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+\n")
            # flip axis
            #Rplotter_file.write("coord_flip()+\n")

            # centre plot tile
            MainRplotter_file.write("theme(plot.title=element_text(hjust=0.5))+\n") 

            # make tick marks longer
            MainRplotter_file.write("theme(axis.ticks.length=unit(1,'cm'))+\n") 

            # add space to top and right hand side of plot
            MainRplotter_file.write("theme(plot.margin=unit(c(t=30,r=30,b=0,l=0),'mm'))+\n")
            

            # add plot title
            # Rplotter_file.write(f"ggtitle('Placeholder')+\n") 
        
            # fix x axis lenght based on scaffold length rather than the end of the data
            MainRplotter_file.write(f"scale_x_continuous(limits=c(0,{ScaffoldMetricsDf[Scaf]['ScaffoldLength']}),labels = label_number(scale= 1e-6))\n") 
            
            # define aspect ratio of plot
            # MainRplotter_file.write("coord_fixed()\n")
            # set limits on x axis
            # Rplotter_file.write(f"scale_x_discrete(limits=c(0, {ScaffoldMetricsDf[Scaf]['ScaffoldLength']}))\n")
    # set a working directory
    MainRplotter_file.write(f"pdf('CombinedPlot_{args.sw}_{args.bp_size}_Win.pdf', width = 100, height = 100)\n" )                   
    MainRplotter_file.write(f'ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8, labels = c("A", "B", "C","D","E","F","G","H"),font.label = list(size = 56, color = "red"),ncol = 1, nrow = 8)\n') 
    MainRplotter_file.write("dev.off()\n")

    # run R script
    # os.system(f"Rscript {os.path.abspath(ScriptsDirPath)}/Scaf_{Scaf.lstrip('RL_')}_Plot.R")

os.system(f"Rscript {os.path.abspath(ScriptsDirPath)}/MainCombinedPlot.R")



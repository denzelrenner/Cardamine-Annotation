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

vars.add_argument('--gff_ms',dest='gff_ms', type=str, required=True, help='a gff file produced by barrnap for main scaffolds (MS)')

vars.add_argument('--gff_db',dest='gff_db', type=str, required=True, help='a gff file produced by barrnap for debris scaffolds (DB)')

# vars.add_argument('-fa','--fasta',dest='fasta', type=str, required=True, help='the fasta file containing all the scaffolds')

vars.add_argument('-d','--directory',dest='dir', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

# vars.add_argument('-sw','--sliding_window',dest='sw', type=int, required=False,default=1000, help='size for the sliding window within which metrics are calculated. (default=1000)')

# vars.add_argument('-bp','--bp_size',dest='bp_size', type=str, required=False,default='b', help='sliding window value will be interpereted as being in b, kb, or mb. (default=b)')

vars.add_argument('--flip',dest='flip_coors', action='store_true',required=False,default=False, help='set this argument (i.e --flip) to flip the coordinates of the dot plot')

vars.add_argument('--ignore_zero',dest='ignore_zero', action='store_true',required=False,default=False, help='set this argument (i.e --ignore_zero) to stop density values of 0 from being plotted')

vars.add_argument('--query_loci',dest='query_loci', nargs='+', type=str, required=False,default='all_loci', help='give the name of the rDNA loci repeats in order of their appearance in the repeat block (i.e 5s or 5s 28s)')

vars.add_argument('--combined_plot',dest='combined_plot', action='store_true',required=False,default=False, help='set this argument to get a combined plot with eahc scaffold')

args = vars.parse_args()

# print(args.query_loci)
# sys.exit()


# sorting out the name of output plots
ChosenLoci = None

if args.query_loci == 'all_loci':
    ChosenLoci = 'all_loci'

else:
    ChosenLoci = '_'.join(args.query_loci)

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

# global variable which has Scaffold metrics. format {'Scaffold':{'length(in bp)':200000}}
ScaffoldMetricsDf = {}

# create data frame which we will add data to
ArryOfOutputDataDicts = []

# store data in pandas df
GffOutputMS = pd.read_csv(filepath_or_buffer=f"{args.gff_ms}",sep='\t',skiprows=1,header=None)
GffOutputDB = pd.read_csv(filepath_or_buffer=f"{args.gff_db}",sep='\t',skiprows=1,header=None)

# for the MainScaffolds 
# grab the loci information from column 8
# for ProductLine in GffOutputMS[8]:


#     pass
# get all the different loci in the scaffold
# if all the expected loci are not in AllLoci varibale, then add them with a value of 0 to the final table
AllLoci = list(set([grep_product(line=ProductLine) for ProductLine in GffOutputMS[8]]))

# loop through all loci and get their quantity in the dataset
for Loci in AllLoci:
# for Loci in args.query_loci:

    # main scaffold variables
    NumberOfPartiallyAlignedMS = 0
    NumberOfFullyAlignedMS = 0
    TotalNumberMS = 0

    # debris variables
    NumberOfPartiallyAlignedDB = 0
    NumberOfFullyAlignedDB = 0
    TotalNumberDB = 0

    for ProductLine in GffOutputMS[8]:

        # try to match the loci in the line
        FoundHit = re.search(rf'{Loci}',ProductLine)

        if FoundHit:

            # smth += 1
            TotalNumberMS += 1
        
            # try to see if it is whole or partially aligned
            IsAligned = re.search(rf'aligned',ProductLine)

            if IsAligned:

                # smth += 1
                NumberOfPartiallyAlignedMS += 1

    # now repeat for the debris as well
    for ProductLine in GffOutputDB[8]:

        # try to match the loci in the line
        FoundHit = re.search(rf'{Loci}',ProductLine)

        if FoundHit:

            # smth += 1
            TotalNumberDB += 1
        
            # try to see if it is whole or partially aligned
            IsAligned = re.search(rf'aligned',ProductLine)

            if IsAligned:

                # smth += 1
                NumberOfPartiallyAlignedDB += 1
    
    # # row for whole loci i.e fully aligned
    # ArryOfOutputDataDicts.append(
    #     {'LociType':f'{Loci}_Whole',
    #      'Debris':TotalNumberDB-NumberOfPartiallyAlignedDB,
    #      'MainScaffold':TotalNumberMS-NumberOfPartiallyAlignedMS,
    #      'Combined':(TotalNumberDB-NumberOfPartiallyAlignedDB)+(TotalNumberMS-NumberOfPartiallyAlignedMS)}
    # )

    # # row for partially aligned loci
    # ArryOfOutputDataDicts.append(
    #     {'LociType':f'{Loci}_Partial',
    #      'Debris':NumberOfPartiallyAlignedDB,
    #      'MainScaffold':NumberOfPartiallyAlignedMS,
    #      'Combined':NumberOfPartiallyAlignedMS+NumberOfPartiallyAlignedDB})
    
    ArryOfOutputDataDicts.append(
        {'ScaffoldType':'MainScaffolds',
         'rDNAType':Loci,
         'Alignment':'Complete',
         'rDNAQuantity':TotalNumberMS-NumberOfPartiallyAlignedMS}
    )

    ArryOfOutputDataDicts.append(
        {'ScaffoldType':'MainScaffolds',
         'rDNAType':Loci,
         'Alignment':'Partial',
         'rDNAQuantity':NumberOfPartiallyAlignedMS}
    )

    ArryOfOutputDataDicts.append(
        {'ScaffoldType':'Debris',
         'rDNAType':Loci,
         'Alignment':'Complete',
         'rDNAQuantity':TotalNumberDB-NumberOfPartiallyAlignedDB}
    )

    ArryOfOutputDataDicts.append(
        {'ScaffoldType':'Debris',
         'rDNAType':Loci,
         'Alignment':'Partial',
         'rDNAQuantity':NumberOfPartiallyAlignedDB}
    )
# ScaffoldType      rDNAType  rDNAQuality  rDNAQuantity
# Main Scaffolds     18S      Whole         30
# Main sacffolds     18S      Align         10
# Main Scaffolds
# create intermediate df
OutputDataFrame = pd.DataFrame(ArryOfOutputDataDicts)

print(OutputDataFrame)
# sys.exit()




# for the output table we are making we want to have every type of loci
# 5S 5.8S 18S 28S
# loci type                  debris mainscaffolds total
#5S Whole
#5S partially aligned
#5.8S whole
#5.8S partially aligned

# save table in tsv
OutputDataFrame.to_csv(path_or_buf=f'{OutputDataDirPath}/Final_Table.tsv',sep='\t',index=False)

# open plotting file
with open(f"{ScriptsDirPath}/MainPlot.R",'w') as MainRplotter_file:

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
    MainRplotter_file.write(f"df <- read.table(file='{os.path.abspath(OutputDataDirPath)}/Final_Table.tsv',header=TRUE,sep='\\t')\n")

    # start graphics driver
    # load graphics driver
    # MainRplotter_file.write(f"pdf('Scaf_{Scaf}_Plot_{args.sw}_{args.bp_size}_Win.pdf', width = 25, height = 25)\n")
    
    # reorder legend
    #convert 'position' to factor and specify level order

    MainRplotter_file.write("df$Alignment <- factor(df$Alignment, levels=c('Partial', 'Complete'))\n")

    # main plot code ScaffoldType      rDNAType  rDNACompleteness  rDNAQuantity
    MainRplotter_file.write(f"p <- ggplot(df) +\n")
    # MainRplotter_file.write("geom_bar(aes(x=rDNAType,y=rDNAQuantity,fill=rDNACompleteness),position='stack',stat='identity') +\n")
    MainRplotter_file.write("geom_bar(aes(x=ScaffoldType,y=rDNAQuantity,fill=Alignment),colour='black',position='stack',stat='identity') +\n")

    MainRplotter_file.write("facet_grid(~rDNAType)+\n")
    MainRplotter_file.write(f"labs(x = 'Scaffold Type',y = 'Number of Loci') +\n")
    
    # # remove grid from background of plot
    # # remove grid from background of plot
    MainRplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
    MainRplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
    MainRplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
    MainRplotter_file.write('''theme(axis.title.y = element_text(vjust=-2.5,margin = unit(c(t=0, r=7, b=0, l=0), "mm")),
        axis.title.x = element_text(vjust=5.5,margin = unit(c(t=10, r=0, b=0, l=0), "mm")),axis.text = element_text(),
        legend.text=element_text(),axis.title=element_text(),
        legend.title = element_text(),legend.key.size = unit(10, 'mm'),
        axis.ticks.y = element_line(), axis.line.y = element_line(),
        axis.ticks.x = element_line(),axis.line.x = element_line())+\n''')
    
    # Use custom colors
    MainRplotter_file.write("scale_fill_manual(values=c('#999999','#E69F00'))\n")
 

    # # centre plot tile
    # MainRplotter_file.write("theme(plot.title=element_text(hjust=0.5))+\n") 

    # # make tick marks longer
    # MainRplotter_file.write("theme(axis.ticks.length=unit(1,'cm'))+\n") 

    # # add space to top and right hand side of plot
    # MainRplotter_file.write("theme(plot.margin=unit(c(t=70,r=50,b=0,l=0),'mm'))+\n")
    # # add plot title
    # # Rplotter_file.write(f"ggtitle('Placeholder')+\n") 
    # # fix x axis lenght based on scaffold length rather than the end of the data
    # MainRplotter_file.write(f"scale_x_continuous(limits=c(0,{ScaffoldMetricsDf[Scaf]['ScaffoldLength']}),labels = label_number(scale= 1e-6))\n") 

    # finish off script for combined plot now that we have added information for each scaffold
    MainRplotter_file.write(f"pdf('FinalPlot.pdf', width = 10, height = 10)\n" )                   
    MainRplotter_file.write(f'ggarrange(p, labels = c("A"),font.label = list(size = 10, color = "red"),ncol = 1, nrow = 1)\n') 
    MainRplotter_file.write("dev.off()\n")

    # # run R script
    # os.system(f"Rscript {os.path.abspath(ScriptsDirPath)}/Scaf_{Scaf.lstrip('RL_')}_Plot.R")

# only run the combined script if specified
os.system(f"Rscript {os.path.abspath(ScriptsDirPath)}/MainPlot.R")



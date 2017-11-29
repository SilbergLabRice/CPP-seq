######################
###CPP-seq Analysis###
######################

########
#Inputs:
########
    #ParentSequences folder = (1) A FASTA file of variable region (Gene+linker) with the gene name as the ID (e.g. TnAK), and (2) FASTA file of Transposon used for CPP-Seq
    #MiSeq_Data folder: For each sequence in the variable region FASTA file you need two FASTQ files (1) from unselected library, and (2) from selected library 

#########
#Outputs:
#########
    #(1) Summary text file for all sequences analyzed, for each seqeuence (2) FoldChange and (3)FishersTest and the following for both the selected and unselected libraries of each sequence: (2)Raw Counts, (3) Tagged counts, (4)NuclotidePositionCounts, (5)AminoAcidTerminiCounts 

import re, sys
import UTIL.UTIL as util
import Bio.SeqIO
import sys


#Summary File output 
print('Processing Sequences . . .')
f = open("Outputs/Batch Process Summary.txt", 'w')
sys.stdout = f

gene_names = []
variable = 'Data/variable.fasta'
transposon = 'Data/transposon.fasta'


for record in Bio.SeqIO.parse(variable, "fasta"):
    gene_names.append(record.id)

kmer_param = 54        #parameter for how long of a transpson tag (bp directly adjacent to VARIABLE region) to search for, requires a length such that they are unique sequences

#Transposon 5' and 3' Tags      
for record in Bio.SeqIO.parse(transposon, "fasta"):
    print('Determining 5 and 3 prime tags from the ' + str(record.id) + ' transposon:')
    
    transposon = record.seq.upper()
    five_fw_tag = transposon[len(transposon)-kmer_param:len(transposon)]
    five_rev_tag =  five_fw_tag.reverse_complement()
    three_fw_tag = transposon[0:kmer_param]
    three_rev_tag = three_fw_tag.reverse_complement()
    
five_fw_tag = str(five_fw_tag)
five_rev_tag = str(five_rev_tag)
three_fw_tag = str(three_fw_tag)
three_rev_tag = str(three_rev_tag)

    
#Process raw sequencing data to extract all gene-transposon junctions and their position in the gene
for gene in gene_names:
    
    gene_name = gene 
    
    #Read sequence from variable.fasta
    for record in Bio.SeqIO.parse(variable, "fasta"):
        if record.id == gene_name:
            gene_fwd = str(record.seq)
            gene_rev = str(record.seq.reverse_complement())
    
    #GenerateLibraryofPossibleInsertionSites
    util.GenerateTransposonRecognitionSites(gene_name, gene_fwd, gene_rev)
    
    #Analyzing unselected and selected counts 
    filenames = ['Data/' + gene_name + '_unselected.fastq', 'Data/' + gene_name + '_selected.fastq']
    
    for file in filenames:
        #RAW data processing
        util.TransposonGeneJunctionExtraction(file, five_fw_tag, five_rev_tag, three_fw_tag, three_rev_tag, kmer=kmer_param)
        util.AdjacentMerExtraction(file, five_fw_tag, five_rev_tag, three_fw_tag, three_rev_tag, kmer=kmer_param)
        util.PositionExtraction(file, gene_name, gene_fwd, gene_rev)
        
        #Transposon Frequency Analysis (TFA)
        util.NucleotideInsertionSiteCounter(file, gene_fwd)
        util.AminoAcidStartSiteCounter(file, gene_fwd)
        
    #Calculate Selected/Unselected Fold ratios
    util.FoldChangeCalculator(gene_name)
    
    #Statistics using Parallel v. Antiparallel Reference
    #Fisher's exact test (requires reads at every positon)
    util.FisherTest(gene_name)
    #Fold Change test (estimates variance from AP reads)
    
f.close() #close summary file 

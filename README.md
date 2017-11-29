CPP-seq

A custom python pipeline for analyzing Illumina MiSeq data for Circular Permutation Profiling. 


Inputs:
    
    Inputs are stored in the 'Data' directory 
    
    Inputs consist of the following:
    (1) A FASTA file of variable region (Gene+linker) with the gene name as the ID (e.g. TnAK)
    (2) FASTA file of Transposon used for CPP-Seq
    
    MiSeq_Data folder: For each sequence in the variable region FASTA file you need two FASTQ files 
    (1) from unselected library
    (2) from selected library 


Outputs:
    
    Outputs are stored in the 'Outpus' directory
    
    Outputs consist of:
    (1) Summary text file for all sequences analyzed, for each seqeuence 
    (2) FoldChange
    (3) FishersTest 
    as well as one copy of the following for both the selected and unselected libraries of each sequence:
    (2) Raw Counts
    (3) Tagged counts
    (4) NuclotidePositionCounts
    (5) AminoAcidTerminiCounts 

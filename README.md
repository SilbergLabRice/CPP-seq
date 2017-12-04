# CPP-seq

A custom python pipeline for analyzing Illumina MiSeq data for Circular Permutation Profiling. 

## Script overview  
### CPP-SequenceProcessing.py 
    Processes the unselected and selected library MiSeq data and counts the insertions at every position within the gene(s) of interest being analayzed.

    Inputs:

        Inputs are stored in the 'Data' directory 

        Inputs consist of the following:
        (1) A FASTA file of variable region (gene+linker, 5'-start codon->linker-3') with the gene name as the ID (e.g. TnAK)
        (2) FASTA file of transposon used for CPP-Seq (e.g. MT2)

        For each sequence in the variable region FASTA file you need two MiSeq FASTQ files 
        (1) from unselected library (e.g. TnAK_unselected.fastq)
        (2) from selected library (e.g. TnAK_selected.fastq)

        Example Variable region and Transposon FASTA files are included in the 'Data' directory 
        Example MiSeq data can be found at: 
            TnAK_unselected.fastq (https://www.ncbi.nlm.nih.gov/sra/SRX3427327)
            TnAK_selected.fastq (https://www.ncbi.nlm.nih.gov/sra/SRR6327684)

    Outputs:

        Outputs are stored in the 'Outputs' directory

        Outputs consist of (11 files/gene):
        For each gene analyzed (3/gene)
        (1) Summary text file for all sequences analyzed, for each seqeuence 
        (2) FoldChange
        (3) FishersTest 

        as well as one copy of the following for both the selected and unselected libraries of each sequence (8 total/gene):

        (1) Raw Counts
        (2) Tagged counts
        (3) NuclotidePositionCounts
        (4) AminoAcidTerminiCounts 
       
### CPP-Plots.py 
    Uses the processed sequencing data to plot insertional frequencies with respect to the protein sequence yielding plots in the 'Plot' directory for each possible reading frame. 

##  Authors

* Joshua T Atkinson - *NGS processing and analysis* 
* Quan Zhou - *Statistical analysis*

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

* Simon Anders and Wolfgang Huber for the DESeq program package (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) that was the basis of our statistical analysis. 


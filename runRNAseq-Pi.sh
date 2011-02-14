#!/bin/bash

# runRNaseq-Pi.sh
# Example of RNAseq-Pi script/modules filter functions, all this functions will
# be controlled in the real pipeline.
# Copyright (C) 2011 Juan Caballero [Institute for Systems Biology]

# add Bowtie in your PATH
export PATH=$PATH:/path/to/bowtie
# add our Perl scripts in your PATH (check permisions first)
export PATH=$PATH:/path/to/scripts
# NPROC -> number of CPUs to use
NPROC=4 
# location of the indexes
BWINDEX=/path/to/bowtie/indexes
# Bowtie aditional parameters
BWPARAM="--quiet -p $NPROC -S --sam-nohead -k 1 -v 2"

for FQ in *.fq # Process all the FASTQ files in the current directory
do
    ID=${FQ%fq}
    # cascade mode, if you need to change step order, must change outputs/inputs
    
    # step1 -> filter low quality reads
    filterQuality.pl -v -i $FQ -o $ID.qual_OK.fq -b $ID.qual_BAD.fq
    
    # step2 -> remove vector matches
    bowtie $BWPARAM $BWINDEX/UniVec_Core -q $ID.qual_OK.fq -un $ID.vector_OK.fq $ID.bowtie_vector_BAD.sam
    
    # step3 -> filter low complexity reads
    filterLowComplex.pl -v -i $ID.vector_OK.fq -o $ID.complex_OK.fq -b $ID.complex_BAD.fq
    
    # step4 -> remove rRNA/MT matches
    bowtie $BWPARAM $BWINDEX/human.GRCh37.61.rRNA-MT -q $ID.complex_OK.fq -un $ID.riboMT_OK.fq $ID.bowtie_riboMT_BAD.sam
    
    # step5 -> remove repeats matches
    bowtie $BWPARAM $BWINDEX/human_RepBase15.10 -q $ID.riboMT_OK.fq -un $ID.repeat_OK.fq $ID.bowtie_repeat_BAD.sam
    
    # step6 -> remove ERCC matches
    bowtie $BWPARAM $BWINDEX/ERCC_reference_081215 -q $ID.repeat_OK.fq -un $ID.ercc_OK.fq $ID.bowtie_ercc_BAD.sam
    
    # collect the filtered reads
    mv $ID.ercc_OK.fq $ID.FILTERED.fq
    
    # remove/compress temporary files
    for FILE in *OK.fq *BAD.fq *BAD.sam
    do
        rm $FILE
        #gzip $FILE
    done
    
    # step7 -> assembly, reference aligment, etc
    
done

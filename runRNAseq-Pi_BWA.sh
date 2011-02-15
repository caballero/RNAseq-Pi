#!/bin/bash

# runRNaseq-Pi_BWA.sh
# Example of RNAseq-Pi script/modules filter functions, all this functions will
# be controlled in the real pipeline.
# Copyright (C) 2011 Juan Caballero [Institute for Systems Biology]

# add BWA to your PATH
export PATH=$PATH:/proj/hoodlab/share/programs/bwa
# add our Perl scripts in your PATH (check permisions first)
export PATH=$PATH:/proj/hoodlab/share/programs/RNAseq-Pi/bin
# NPROC -> number of CPUs to use
NPROC=6
# location of the indexes
BWINDEX=/proj/hoodlab/share/programs/RNAseq-Pi/data/bwa-indexes
# BWA aditional parameters
BWPARAM="-t $NPROC"

for FQ in *.fq # Process all the FASTQ files in the current directory
do
    echo "processing $FQ"
    ID=${FQ%.fq}
    # cascade mode, if you need to change step order, must change outputs/inputs

    # step1 -> filter low quality reads
    echo "  #1 low quality filter" 
    date
    filterQuality.pl -v -i $FQ -o $ID.qual_OK.fq -b $ID.bwa_qual_BAD.fq
    date

    # step2 -> remove vector matches
    echo "  #2 vector filter"
    date    
    bwa aln $BWPARAM $BWINDEX/solexa_primers $ID.qual_OK.fq > $ID.bwa_vector.sai
    bwa samse $BWINDEX/solexa_primers $ID.qual_OK.fq $ID.bwa_vector.sai > $ID.bwa_vector.sam
    getUnmapSam.pl -i $ID.bwa_vector.sam -o $ID.bwa_vector_OK.fq -m $ID.bwa_vector_BAD.sam
    rm $ID.bwa_vector.sam $ID.bwa_vector.sai
    date

    # step3 -> filter low complexity reads
    echo "  #3 low complex filter"
    date
    filterLowComplex.pl -v -i $ID.vector_OK.fq -o $ID.complex_OK.fq -b $ID.bwa_complex_BAD.fq
    date

    # step4 -> remove rRNA/MT matches
    echo "  #4 rRNA+MT filter"    
    date
    bwa aln $BWPARAM $BWINDEX/human.GRCh37.61.rRNA-MT $ID.complex_OK.fq > $ID.bwa_riboMT.sai
    bwa samse $BWINDEX/human.GRCh37.61.rRNA-MT $ID.complex_OK.fq $ID.bwa_riboMT.sai > $ID.bwa_riboMT.sam
    getUnmapSam.pl -i $ID.bwa_riboMT.sam -o $ID.bwa_riboMT_OK.fq -m $ID.bwa_riboMT_BAD.sam
    rm $ID.bwa_riboMT.sam $ID.bwa_riboMT.sai 
    date

    # step5 -> remove repeats matches
    echo "  #5 repeats filter"    
    date
    bwa aln $BWPARAM $BWINDEX/human_RepBase15.10 $ID.riboMT_OK.fq > $ID.bwa_repeat.sai
    bwa samse $BWINDEX/human_RepBase15.10 $ID.riboMT_OK.fq $ID.bwa_repeat.sai > $ID.bwa_repeat.sam
    getUnmapSam.pl -i $ID.bwa_repeat.sam -o $ID.bwa_repeat_OK.fq -m $ID.bwa_repeat_BAD.sam
    rm $ID.bwa_repeat.sam $ID.bwa_repeat.sai 
    date

    # step6 -> remove ERCC matches
    echo "  #6 ERCC filter"    
    date
    bwa aln $BWPARAM $BWINDEX/ERCC_reference_081215 $ID.repeat_OK.fq > $ID.bwa_ercc.sai
    bwa samse $BWINDEX/ERCC_reference_081215 $ID.repeat_OK.fq $ID.bwa_ercc.sai > $ID.bwa_ercc.sam
    getUnmapSam.pl -i $ID.bwa_ercc.sam -o $ID.bwa_ercc_OK.fq -m $ID.bwa_ercc_BAD.sam
    rm $ID.bwa_ercc.sam $ID.bwa_ercc.sai 
    date

    # collect the filtered reads
    mv $ID.ercc_OK.fq $ID.bwa_FILTERED.fq
    
    # remove/compress temporary files
    echo "  compressing files"
    date
    for FILE in *_BAD.fq *_BAD.sam $ID.bwa_FILTERED.fq
    do
        gzip $FILE
    done
    date

    echo "  removing temp files"
    rm *_OK.fq

    # step7 -> assembly, reference aligment, etc

done


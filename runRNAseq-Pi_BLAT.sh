#!/bin/bash

# runRNaseq-Pi_BLAT.sh
# Example of RNAseq-Pi script/modules filter functions, all this functions will
# be controlled in the real pipeline.
# Copyright (C) 2011 Juan Caballero [Institute for Systems Biology]

# add Blat to your PATH
export PATH=$PATH:/proj/hoodlab/share/programs/blat
# add our Perl scripts in your PATH (check permisions first)
export PATH=$PATH:/proj/hoodlab/share/programs/RNAseq-Pi/bin
# NPROC -> number of CPUs to use
NPROC=6
# location of the indexes
BWINDEX=/proj/hoodlab/share/programs/RNAseq-Pi/data/blat-indexes
# BWA aditional parameters
BWPARAM="-fastMap -q=DNA -t=DNA"

for FQ in *.fq # Process all the FASTQ files in the current directory
do
    echo "processing $FQ"
    ID=${FQ%.fq}
    # cascade mode, if you need to change step order, must change outputs/inputs

    # step0 -> convert FastQ to Fasta
    echo "  #0 converting fastq->fasta"
    date
    cat $FQ | fq2fa.pl > $ID.fa
    date
    
    # step1 -> filter low quality reads
    echo "  #1 low quality filter" 
    date
    filterQuality.pl -v -f fa -i $ID.fa -o $ID.qual_OK.fa -b $ID.blat_qual_BAD.fa
    date

    # step2 -> remove vector matches
    echo "  #2 vector filter"
    date   
    blat $BWINDEX/solexa_primers.2bit $ID.qual_OK.fa $BWPARAM $ID.blat_vector_BAD.psl
    removeBlatHit.pl $ID.blat_vector_BAD.psl $ID.qual_OK.fa > $ID.vector_OK.fa
    date

    # step3 -> filter low complexity reads
    echo "  #3 low complex filter"
    date
    filterLowComplex.pl -v -f fa -i $ID.vector_OK.fa -o $ID.complex_OK.fa -b $ID.blat_complex_BAD.fa
    date

    # step4 -> remove rRNA/MT matches
    echo "  #4 rRNA+MT filter"    
    date
    blat $BWINDEX/human.GRCh37.61.rRNA-MT.2bit $ID.complex_OK.fa $BWPARAM $ID.blat_riboMT_BAD.psl
    removeBlatHit.pl $ID.blat_riboMT_BAD.psl $ID.complex_OK.fa > $ID.riboMT_OK.fa
    date

    # step5 -> remove repeats matches
    echo "  #5 repeats filter"    
    date
    blat $BWINDEX/human_RepBase15.10.2bit $ID.riboMT_OK.fa $BWPARAM $ID.blat_repeat_BAD.psl
    removeBlatHit.pl $ID.blat_repeat_BAD.psl $ID.riboMT_OK.fa > $ID.repeat_OK.fa
    date

    # step6 -> remove ERCC matches
    echo "  #6 ERCC filter"    
    date
    blat $BWINDEX/ERCC_reference_081215.2bit $ID.repeat_OK.fa $BWPARAM $ID.blat_ercc_BAD.psl
    removeBlatHit.pl $ID.blat_ercc_BAD.psl $ID.repeat_OK.fa > $ID.ercc_OK.fa
    date

    # collect the filtered reads
    mv $ID.ercc_OK.fa $ID.blat_FILTERED.fa
    
    # remove/compress temporary files
    echo "  compressing files"
    date
    for FILE in *_BAD.fa *_BAD.psl $ID.blat_FILTERED.fa
    do
        gzip $FILE
    done
    date

    echo "  removing temp files"
    rm *_OK.fa

    # step7 -> assembly, reference aligment, etc

done


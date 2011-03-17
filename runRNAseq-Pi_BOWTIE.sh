#!/bin/bash

# runRNaseq-Pi.sh
# Example of RNAseq-Pi script/modules filter functions, all this functions will
# be controlled in the real pipeline.
# Copyright (C) 2011 Juan Caballero [Institute for Systems Biology]

# add our Perl scripts in your PATH (check read/exec permisions first)
export PATH=$PATH:/proj/hoodlab/share/programs/RNAseq-Pi/bin

# NPROC -> number of CPUs to use
NPROC=4

# Compressor to use
ZIP=gzip

# Location of the indexes
BWINDEX=/proj/hoodlab/share/programs/indexes

# Bowtie
export PATH=$PATH:/path/to/bowtie
BWPARAM="--quiet -p $NPROC -S --sam-nohead -k 1 -v 2 -q"

# MapSplice
MAPSPLICE=/proj/hoodlab/share/programs/mapsplice/bin/mapsplice_segments.py
GENOME=hs37.61
GENOMEDIR=/proj/hoodlab/share/programs/Ensembl
READSIZE=75
MSPARAM="-c $GENOMEDIR/$GENOME.fa -B $BWINDEX/$GENOME -t $GENOMEDIR/$GENOME.pseudo.gtf -w $READSIZE -L 25 -Q fq -X $NPROC --full-running --semi-canonical"

# Cufflinks
export PATH=$PATH:/proj/hoodlab/share/programs/cufflinks
CFPARAM="-r $GENOMEDIR/$GENOME.gtf -s $GENOMEDIR/$GENOME.fa"

for FQ in *.fq # Process all the FASTQ files in the current directory
do
    ID=${FQ%fq}
    # cascade mode, if you need to change step order, must change outputs/inputs
    
    # step1 -> filter low quality reads
    filterQuality.pl -i $FQ -o $ID.qual_OK.fq -b $ID.bowtie_qual_BAD.fq
    
    # step2 -> remove vector matches
    bowtie $BWPARAM --un $ID.vector_OK.fq $BWINDEX/UniVec_Core $ID.qual_OK.fq | perl -lane 'print unless($F[1] == 4)' > $ID.bowtie_vector_BAD.sam
    
    # step3 -> filter low complexity reads
    filterLowComplex.pl -i $ID.vector_OK.fq -o $ID.complex_OK.fq -b $ID.bowtie_complex_BAD.fq
    
    # step4 -> remove rRNA/MT matches
    bowtie $BWPARAM --un $ID.riboMT_OK.fq $BWINDEX/human.GRCh37.61.rRNA-MT $ID.complex_OK.fq | perl -lane 'print unless($F[1] == 4)' > $ID.bowtie_riboMT_BAD.sam
    
    # step5 -> remove repeats matches
    bowtie $BWPARAM --un $ID.repeat_OK.fq $BWINDEX/human_RepBase15.10 $ID.riboMT_OK.fq | perl -lane 'print unless($F[1] == 4)' > $ID.bowtie_repeat_BAD.sam
    
    # step6 -> remove ERCC matches
    bowtie $BWPARAM --un $ID.ercc_OK.fq $BWINDEX/ERCC_reference_081215 $ID.repeat_OK.fq | perl -lane 'print unless($F[1] == 4)' > $ID.bowtie_ercc_BAD.sam
    
    # collect the filtered reads
    mv $ID.ercc_OK.fq $ID.bowtie_FILTERED.fq
    
    # compress temporary files
    for FILE in *BAD.fq *BAD.sam
    do
        $ZIP $FILE
    done
    
    # remove temporary files
    rm *_OK.fq
    
    # step7 -> alignment with MapSplice
    python $MAPSPLICE $MSPARAM -u $ID.bowtie_FILTERED.fq -o $ID.mapsplice
    grep chr $ID.mapsplice/alignments.sam | perl -lane '$F[1] == 16 ? print "$_\tXS:A:-" : print "$_\tXS:A:+"' > $ID.mapsplice.sam
    # if we ran using fusion detection, uncomment the next line to keep the file
    #mv $ID.mapsplice/fusion_junction $ID.mapsplice.fusion_junction
    rm -rf $ID.mapsplice
    
    # step8 -> Cufflink assembly and transcript quantification
    sort -k 3,3 -k 4,4n $ID.mapsplice.sam > $ID.mapsplice.sorted.sam
    rm $ID.mapsplice.sam
    mv $ID.mapsplice.sorted.sam $ID.mapsplice.sam
    cufflinks -q -p $NPROC -o $ID.cufflinks $ID.mapsplice.sam
    cuffcompare $CFPARAM -o $ID.cuff $ID.cufflinks/transcripts.gtf
    # transcript annotation and quantification is in $ID.cufflinks/transcripts.tmap
    
    $ZIP $ID.mapsplice.sam
    
    # END
done

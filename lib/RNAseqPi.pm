package RNAseqPi;

=head1 NAME

RNAseqPi

=head1 DESCRIPTION

Perl module for RNAseq-Pi.

=head1 AUTHOR

Juan Caballero
Institute for Systems Biology @ 2010

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
our $VERSION = '1.00';

=head2 

Sub: loadParam

Description: Load a configuration file

Call: loadParam($conif)

Return: 1 => ok

=cut

sub loadParam {
    my ($config) = @_;
    open C, "$config" or die "cannot open $config\n";
    while (<C>) {
        chomp;
        next if (m/^#/);
        my ($name, $value) = split (/\t/, $_);
        next unless (defined $name and defined $value);
        if ($name eq 'step') {
            push @{ $param{'step'} }, $value;
        }
        else {
            $param{$name} = $value;
        }
    }
    close C;
    
    return 1;
}

=head2 

Sub: writeFastq

Description: write a Fastq file 

Call: writeFastq($ref_seq, $fastq);

Return: 1 => OK

=cut

sub writeFastq {
    my ($ref_seq, $fastq) = @_;
    open FQ, ">$fastq" or die "cannot write $fastq\n";
    while (my $seqid = each %{ $ref_seq }) {
        my $seq  = $str_seq{$seqid}{'seq'};
        my $qual = $str_seq{$seqid}{'qual'};
        print FQ "\@$seqid\n$seq\n+\n$qual\n";
    }
    close FQ;
    return 1;
}

=head2 

Sub: writeFasta

Description: write a Fasta file 

Call: writeFasta($ref_seq, $fasta);

Return: 1 => OK

=cut

sub writeFasta {
    my ($ref_seq, $fasta) = @_;
    open FA, ">$fasta" or die "cannot write $fasta\n";
    while (my $seqid = each %{ $ref_seq }) {
        my $seq  = $str_seq{$seqid}{'seq'};
        print FA ">$seqid\n$seq\n";
    }
    close FA;
    return 1;
}

=head2 

Sub: mapSequences

Description: run the short-read mapper (BWA, Bowtie or Blat, in future any SAM compatible mapper)

Call: mapSequences($mapper, $reference, $sam, $key);

Return: 1 => OK

=cut

sub mapSequences {
    my ($mapper, $reference, $ref_seq, $sam, $key) = @_;
    
    my $key ||= int(rand 100000);
    
    checkIndex($mapper, $reference);
    
    if ($mapper eq 'bwa') {
        my $fq  = "$key.fq";
        writeFastq($ref_seq, $fq);
        runBWA($reference, $fq, $sam, $map_opt);
        unlink($fq);        
    } elsif ($mapper eq 'bowtie') {
        my $fq  = "seq$key.fq";
        writeFastq($dbname, $fq, $filter);
        runBowtie($reference, $fq, $sam, $map_opt);
        unlink($fq);        
    } elsif ($mapper eq 'blat') {
        my $fa  = "seq$key.fa";
        writeFasta($dbname, $fa, $filter);
        runBlat($reference, $fa, $sam, $map_opt);
        unlink($fa);
    } else {
        die "I don't recognize this mapper: $mapper\n";
    }
    
    return 1;
}

=head2 

Sub: checkIndex

Description: check for index files for a reference sequence

Call: checkIndex($mapper, $reference)

Return: 1 => ok

=cut

sub checkIndex {
    my ($mapper, $reference) = @_;
    my $path_to_idx = $param{'path_to_seqs'}{$mapper};
    my $ok = 1;
    if ($mapper eq 'bowtie') {
        my @suffix = qw(1.ebwt 2.ebwt 3.ebwt 4.ebwt rev.1.ebwt rev.2.ebwt);
        foreach my $sufx (@suffix) {
            $ok = 0 unless (-e "$path_to_idx/$reference.$sufx");
        }
    }
    
    elsif ($mapper eq 'bwa') {
        my @suffix = qw(amb ann bwt pac rbwt rpac rsa sa);
        foreach my $sufx (@suffix) {
            $ok = 0 unless (-e "$path_to_idx/$reference.$sufx");
        }
    }
    
    elsif ($mapper eq 'blat') {
        $ok = 0;
        if (-e "$path_to_idx/$reference.2bit") {
            $ok = 1;
            $param{'mapper_opts'}{'blat_db'} = "$path_to_idx/$reference.2bit";
        } 
        elsif (-e "$path_to_idx/$reference.nib") {
            $ok = 1;
            $param{'mapper_opts'}{'blat_db'} = "$path_to_idx/$reference.nib";
        } 
        else {
            my $fasta_path = $param{'path_to_seqs'}{'fasta'};
            if (-e "$fasta_path/$reference.fa") {
                $ok = 1;
                $param{'mapper_opts'}{'blat_db'} = "$fasta_path/$reference.fa";
            }
        }
    }
    
    else {
        die "I don't recognize this mapper: $mapper\n";
    }
        
    die "indexes for $reference with $mapper are missing\n" if ($ok =! 1);
    
    return 1;
}

=head2 


=head2 

Sub: runBWA

Description: run BWA

Call: runBWA(reference, fasta/fastq, sam_output, bwa_options)

Return: 1 => ok;

=cut

sub runBWA {
    my ($ref, $fq, $sam, $opt) = @_;
    checkIndex('bwa', $ref);
    $opt ||= $param{'mapper_opts'}{'bwa'};
    my $bwa_bin = $param{'path_to_exec'}{'bwa'};
    system("$bwa_bin aln $opt $ref $fq > $fq.$ref.bwa.sai");
    system("$bwa_bin samse $ref $fq.$ref.bwa.sai $fq > $sam");
    unlink("$ref $fq.$ref.bwa.sai");
    return 1;
}

=head2 

Sub: runBowtie

Description: run Bowtie

Call: runBowtie( reference, fasta/fastq, sam_output, bowtie_options)

Return: 1 => ok;

=cut

sub runBowtie {
    my ($ref, $fq, $sam, $opt) = @_;
    checkIndex('bowtie', $ref);
    $opt ||= $param{'mapper_opts'}{'bowtie'};
    my $bowtie_bin = $param{'path_to_exec'}{'bowtie'};
    system("$bowtie_bin $opt $ref $fq > $sam");
    return 1;
}

=head2 

Sub: runBlat

Description: run Blat

Call: runBlat( reference, fasta/fastq, sam_output, blat_options)

Return: 1 => ok;

=cut

sub runBlat {
    my ($ref, $fa, $sam, $opt) = @_;
    checkIndex('blat', $ref);
    $opt ||= $param{'mapper_opts'}{'blat'};
    my $blat_bin = $param{'path_to_exec'}{'blat'};
    my $blat_db  = $param{'mapper_opts'}{'blat_db'};
    my $psl2sam  = $param{'scripts'}{'psl2sam'};
    system("$blat_bin $opt $blat_db $fa > $fa.$ref.blat.psl");
    system("$psl2sam $fa.$ref.blat.psl $sam");
    unlink("$fa.$ref.blat.psl");
    
    return 1;
}

1;

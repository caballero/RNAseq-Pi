#!/usr/bin/perl

=head1 NAME

gtf2fa.pl

=head1 DESCRIPTION

Extract the sequences from a Fasta file using a GTF annotation.

=head1 USAGE

perl gtf2fa.pl <FASTA> <GTF> <OUT>

=head1 EXAMPLES

perl gtf2fa.pl genome.fa genes.gtf genes.fa

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

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

$ARGV[2] or die "usage perl gtf2fa.pl <FASTA> <GTF> <OUT>\n";

my ($ref, $gtf, $out) = @ARGV;
my %chr;
my %seq;
my ($id, $chr, $ini, $end, $dir, $seq, $gen, $trn, $pep);

my %genetic_code = (   
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'X',    # Stop
    'TAG' => 'X',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'X',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
);

warn "reading sequences from $ref\n";
open FA, "$ref" or die "cannot open file $ref\n";
while (<FA>) {
    chomp;
    if (m/^>(.+)/) {
        $id = $1;
        warn "    $id\n";
    }
    else {
        $chr{$id} .= $_;
    }
}
close FA;

warn "reading GTF from $gtf\n";
open GTF, "$gtf" or die "cannot open file $gtf\n";
while (<GTF>) {
    chomp;
    my @a = split (/\t/, $_);
    next unless ($a[2] eq 'exon');
    $chr  = $a[0];
    $ini  = $a[3];
    $end  = $a[4];
    $dir  = $a[6];
    next unless (defined $chr{$chr});
    
    s/CUFF/FEASTSEQ/g;
    m/gene_id "(.+?)"; transcript_id "(.+?)"/;
    $gen = $1;
    $trn = $2;
    #next if ($gen =~ m/^ENS/ and $trn =~ m/^ENS/);
    $seq  = substr ($chr{$chr}, $ini - 1, $end - $ini);
    if ($dir eq '-') {
        $seq = revcomp($seq);
        if (defined $seq{"$gen|$trn"}) {
            $seq{"$gen|$trn"} = $seq . $seq{"$gen|$trn"};
        }
        else {
            $seq{"$gen|$trn"} = $seq;
        }
    }
    else {
        $seq{"$gen|$trn"} .= $seq;
    }
    
}
close GTF;

warn "writing sequences in $out\n";
open O, ">$out" or die "cannot write file $out\n";
while (($id, $seq) = each %seq) {
    $pep = getORF($seq);
    next unless ( (length $pep) > 30);
    print O ">$id\n";
    while ($pep) {
        my $s = substr ($pep, 0, 80);
        print O "$s\n";
        substr ($pep, 0, 80) = '';
    }
}
close O;

sub revcomp {
    my ($s) = @_;
    my  $r  =  reverse $s;
        $r  =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}

sub getORF {
    my ($seq) = @_;
    # FRAME +1
    my $orf1  = translate($seq); 
    # FRAME +2
    substr ($seq, 0, 1) = '';
    my $orf2  = translate($seq);
    # FRAME +3
    substr ($seq, 0, 1) = '';
    my $orf3  = translate($seq);
    
    my $pep   = bestORF($orf1, $orf2, $orf3);
    return $pep;
}

sub translate {
    my ($seq) = @_;
    my  $orf  = undef;
    for(my $i = 0; $i < (length($seq) - 2) ; $i += 3) {
        $orf .= codon2aa(substr($seq, $i, 3));
    }
    my @orf  = split (/X/, $orf);
    my $best = '';
    foreach my $pep (@orf) {
        $best = $pep if ( (length $pep) > (length $best) );
    }
    return $best;
}

sub codon2aa {
    my($codon) = @_;
    $codon = uc $codon;
    if (defined $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }
    else {
        warn "Bad codon \"$codon\"!!\n";
        return '@';
    }
}

sub bestORF {
    my $best = '';
    foreach my $pep (@_) {
        $best = $pep if ( (length $pep) > (length $best) );
    }
    return $best;
}

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
my ($id, $chr, $ini, $end, $dir, $seq, $gen, $trn);

warn "reading sequences from $ref\n";
open FA, "$ref" or die "cannot open file $ref\n";
while (<FA>) {
    chomp;
    (m/^>(.+)/) ? $id = $1 : $chr{$id} .= $_;
}
close FA;

warn "reading GTF from $gtf\n";
open GTF, "$gtf" or die "cannot open file $gtf\n";
while (<GTF>) {
    chomp;
    my @a = split (/\t/, $_);
    next if  ($a[2] ne 'exon');
    $chr  = $a[0];
    $ini  = $a[3];
    $end  = $a[4];
    $dir  = $a[6];
    next if !(defined $chr{$chr});
    
    s/CUFF/FEASTSEQ/g;
    m/gene_id "(.+?)"; transcript_id "(.+?)"/;
    $gen = $1;
    $trn = $2;
    
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
    print ">$id\n";
    while ($seq) {
        my $s = substr ($seq, 0, 80);
        print "$s\n";
        substr ($seq, 0, 80) = '';
    }
}
close O;

sub revcomp {
    my ($s) = @_;
    my  $r  =  reverse $s;
        $r  =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}

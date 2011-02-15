#!/usr/bin/perl -w
use strict;

=head1 NAME

profileSeqs.pl

=head1 DESCRIPTION

Check for bias in base call. This script reads a fastq stream and collect all
nucleotides by positions, the output is a table showing:

POSITION Num_A Num_C Num_G Num_T Num_N

=head1 USAGE

profileSeqs.pl < FASTQ > TABLE

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

my %pos = ();
my $max = 0;
my @nuc = qw/A C G T N/;
my $nuc = join "\t", @nuc;

$/ = "\n\@";
while (<>) {
    s/^\@//;
    my ($id, $seq, $sep, $qual) = split (/\n/, $_);
    my $len = length $seq;
    my @seq = split (//, $seq);
    $max = $len if ($max < $len);
    for (my $i = 0; $i <= $len - 1; $i++) {
        $pos{$i}{$seq[$i]}++;
    }
}

print "Pos\t$nuc\n";
for (my $i = 0; $i <= $max; $i++) {
    print "$i";
    foreach my $b (@nuc) {
        my $c = 0;
        $c = $pos{$i}{$b} if (defined $pos{$i}{$b});
        print "\t$c";
    }
    print "\n";
}

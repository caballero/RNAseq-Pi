#!/usr/bin/perl -w
use strict;

=head1 NAME

profileQuals.pl

=head1 DESCRIPTION

Check for bias in quality scores. This script reads a fastq stream and collect all
scores by positions, the output is a table showing:

POSITION AVE_SCORE STD_SCORE

=head1 USAGE

profileQuals.pl < FASTQ > TABLE

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

$/ = "\n\@";
while (<>) {
    s/^\@//;
    my ($id, $seq, $sep, $qual) = split (/\n/, $_);
    my $len = length $qual;
    my @qual = qual2num($qual);
    $max = $len if ($max < $len);
    for (my $i = 0; $i <= $len - 1; $i++) {
        push @{ $pos{$i} }, $qual[$i];
    }
}

print "POS\tMEAN_SCORE\tSTDEV_SCORE\n";
for (my $i = 0; $i <= $max - 1; $i++) {
    my $mean = mean(@{ $pos{$i} });
    my $std  = stdev($mean, @{ $pos{$i} });
    print "$i\t$mean\t$std\n";
}

# Subroutines

sub qual2num {
    my $q = shift @_;
    my @q = split (//, $q);
    my @s = ();
    foreach my $x (@q) { push @s, ord($x) - 64; }
    return @s;
}

sub mean {
    my $sum = 0;
    my $num = 0;
    foreach my $x (@_) {
        $sum += $x;
        $num++;
    }
    if ($num > 0) {
        return $sum / $num;
    } else {
        return 'NA';
    }
}

sub stdev {
    my $mean = shift @_;
    my $sum  = 0;
    my $num  = 0;
    foreach my $x (@_) {
        $num++;
        my $diff = $mean - $x;
        my $diff2 = $diff * $diff;
        $sum += $diff2;
    }
    if ($num > 1) {
        return sqrt( $sum / ($num - 1) );
    } else {
        return 'NA';
    }
}

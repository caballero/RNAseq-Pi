#!/usr/bin/perl

=head1 NAME

illumina1.8To1.5.pl

=head1 DESCRIPTION

Convert a fastq stream from Illumina Casava 1.8 back to 1.5 version. 
Changes:
1. Read name is reformatted, no spaces, added /1 and /2 flags
2. Quality scores are re-scaled from Phred+33 tp Phred+64 (B is the lower 
   value).

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2011

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

my $i = 0;
while (<>) {
    chomp;
    $i++;
    if ($i == 1) {
        my $n = convertName($_);
        $_ = $n;
    }
    elsif ($i == 4) {
        my $q = convertQual($_);
        $_ = $q;
        $i = 0;
    }
    print "$_\n";
}

sub convertName {
    my $orig = shift @_;
    my ($base, $tag) = split (/\s+/, $orig);
    my @tag  = split (/:/, $tag);
    my $code = $tag[-1];
    my $pair = $tag[0];
    $base .= "#$code/$pair";
    return $base;
}

sub convertQual {
    my $qual = shift @_;
    my @qual = split (//, $qual);
    for (my $x = 0; $x <= $#qual; $x++) {
         my $q = chr( ord($qual[$x]) + 31);
         $q = 'B' if ($q eq '@' or $q eq 'A');
         $qual[$x] = $q;
    }
    return join "", @qual;
}


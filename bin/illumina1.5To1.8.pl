#!/usr/bin/perl

=head1 NAME

illumina1.8To1.5.pl

=head1 DESCRIPTION

Convert a fastq stream from Illumina Casava 1.5 to 1.8 version. 
Changes:
1. Read name is reformatted
2. Quality scores are re-scaled from Phred+64 to Phred+33 (# is the lower 
   value).

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

my $i = 0;
while (<>) {
    chomp;
    $i++;
    if ($i == 1) {
        my $n = convertName($_);
        $_ = $n;
    }
    elsif ($i == 3) {
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
    # ORIGINAL NAME: @HWI-ST185_0360:1:1101:1512:2137#0/2
    # MODIFIED NAME: @HWI-ST185_0360:0:0:1:1101:1512:2137 2:N:0:0
    
    my $orig  = shift @_;
    my $new   = undef;
    my @name  = split (/[:#\/]/, $orig);
    $name[0] .= ':0:0'; # RunID and FlowCellID are unknown
    my $pair  = pop @name;
    my $barc  = pop @name;
    $new      = join ":", @name;
    $new     .= " $pair:N:0:$barc"; # NoFiltered and NoControlBits
    return $new;
}

sub convertQual {
    # Scale is moved from Phred+64 to Phred+33 (# is min)
    my $qual = shift @_;
    my @qual = split (//, $qual);
    for (my $x = 0; $x <= $#qual; $x++) {
         my $q = chr( ord($qual[$x]) - 31);
         $q = '#' if ($q eq '!' or $q eq '"');
         $qual[$x] = $q;
    }
    return join "", @qual;
}


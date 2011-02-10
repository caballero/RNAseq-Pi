#!/usr/bin/perl -w
use strict;

=head1 NAME

testComplex.pl

=head1 DESCRIPTION

=cut

use SeqComplex;

print "#seq\tat\tats\tce\tcl1\tcl2\tcl3\tcl4\tcl5\tcl6\tcm1\tcm2\tcm3\tcm4\tcm5\tcm6\tcpg\tct1\tct2\tct3\tct4\tct5\tct6\tcwf\tcz\tgc\tgcs\tket\tpur\n";

while (<>) {
    chomp;
    my $seq = $_;
    my $len = length $_;
    my %res = runAllMethods($seq, $len);
    my @met = sort(keys %res);
    print "$seq";
    foreach my $m (@met) {
        my $res = shift @{ $res{$m} };
        $res = sprintf("%.4f", $res);
        print "\t$res";
    }
    print "\n";
}

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


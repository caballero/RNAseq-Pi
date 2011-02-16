#!/usr/bin/perl

=head1 NAME

removeBlatHit.pl

=head1 DESCRIPTION

Filter a fasta file from hits in a blat search (PSL format)

=head1 USAGE

removeBlatHit.pl PSL FASTA > FASTA_filtered

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

$ARGV[1] or die "usage: removeBlatHit.pl PSL FASTA > FASTA_filtered\n";
my $psl = shift @ARGV;
my $fa  = shift @ARGV;

my %skip = ();
open PSL, "$psl" or die "cannot open $psl\n";
while (<PSL>) {
    next unless (m/^\d/);
    my @line = split (/\t/, $_);
    my $id   = $line[9];
    $skip{$id} = 1;
}
close PSL;

$/ = "\n>";
open FA, "$fa" or die "cannot open $fa\n";
while (<FA>) {
    s/>//g;
    my ($id, $seq) = split (/\n/, $_);
    next if (defined $skip{$id});
    print ">$id\n$seq\n";
}
close FA;

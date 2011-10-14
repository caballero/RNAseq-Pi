#!/usr/bin/perl

=head1 NAME

filterPseudogeneSequences.pl

=head1 DESCRIPTION

Read a FASTA stream and remove pseudogene sequences (marked in the fasta name).

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

$/ = "\n>";
while (<>) {
    s/>//;
    next if (m/pseudogene/);
    print ">$_";
}

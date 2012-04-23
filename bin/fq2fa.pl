#!/usr/bin/perl -w
use strict;

=head1 NAME

fq2fa.pl

=head1 DESCRIPTION

Convert a Fastq file into a regular Fasta file (losing base quality scores).

=cut

my $ln = 0; # Line number
while (<>) {
    $ln++;
    if ($ln == 1) {
        s/^\@/>/;
        s/\s+/_/;   
    }
    elsif ($ln == 2) {
        # do nothing, its the sequence   
    }
    elsif ($ln == 3) {
        next;
    }
    elsif ($ln == 4) {
        $ln = 0;
        next;
    }
    print $_;
    
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


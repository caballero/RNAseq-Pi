#!/usr/bin/perl

=head1 NAME

=head1 DESCRIPTION

checkPairs.pl

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

$ARGV[3] or die "usage: checkPairs.pl PAIR1.fq PAIR2.fq PAIR1_OUT.fq PAIR2_OUT.fq\n";

my $par1_f = shift @ARGV; 
my $par2_f = shift @ARGV;
my $out1_f = shift @ARGV;
my $out2_f = shift @ARGV;

my %par1 = ();

# File opening
my $par1_h = $par1_f;
$par1_h = "gunzip -c  $par1_f | " if ($par1_f =~ m/gz$/);
$par1_h = "bunzip2 -c $par1_f | " if ($par1_f =~ m/bz2$/);

my $par2_h = $par2_f;
$par2_h = "gunzip -c  $par2_f | " if ($par2_f =~ m/gz$/);
$par2_h = "bunzip2 -c $par2_f | " if ($par2_f =~ m/bz2$/);

open P1,  "$par1_h" or die "cannot open $par1_f\n";
open P2,  "$par2_h" or die "cannot open $par2_f\n";
open O1, ">$out1_f" or die "cannot open $out1_f\n";
open O2, ">$out2_f" or die "cannot open $out2_f\n";

print "loading first pair $par1_f\n";
my $np1 = 0;
my $np2 = 0;
my $npp = 0;
my $sid = undef;
$/ = "\n\@";
while (<P1>) {
    s/\@//g;
    if (m/(.+:\d+:\d+:\d+:\d+)#/) {
        $par1{$1} = $_;
        $np1++;
    }
    else {
        die "error parsing $par1_f -> $_\n";
    }
}

print "$np1 sequences loaded\npairing with second pair $par2_f\n";
while (<P2>) {
    s/\@//g;
    if (m/(.+:\d+:\d+:\d+:\d+)#/) {
        $sid = $1;
        $np2++;
        if (defined $par1{$sid}) {
            $npp++;
            print O1 "\@$par1{$sid}";
            print O2 "\@$_";
        }
    }
    else {
        die "error parsing $par2_f -> $_\n";
    }
}
print "$np2 sequences read\n$npp are valid pairs\nDone\n";

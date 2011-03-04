#!/usr/bin/perl

=head1 NAME

mixKnownTranscripts.pl

=head1 DESCRIPTION

=head1 USAGE


=head1 EXAMPLES


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
use Getopt::Long;
use Pod::Usage;

my $pair1  = undef;
my $pair2  = undef;
my $fusion = undef;
my $output = undef;
my $unpair = undef;

my %reads2 = ();

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    '1|pair1:s'        => \$pair1,
    '2|pair2:s'        => \$pair2,
    'f|fusion:s'       => \$fusion,
    'u|unpair:s'       => \$unpair,
    'o|output:s'       => \$output
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $pair1 and defined $pair2);

# opening files
my $pair1_h = $pair1;
$pair1_h = "gunzip  -c $pair1 | " if ($pair1 =~ m/\.gz$/);
$pair1_h = "bunzip2 -c $pair1 | " if ($pair1 =~ m/\.bz2$/);
open P1, "$pair1_h" or die "cannot open $pair1\n";

my $pair2_h = $pair2;
$pair2_h = "gunzip  -c $pair2 | " if ($pair2 =~ m/\.gz$/);
$pair2_h = "bunzip2 -c $pair2 | " if ($pair2 =~ m/\.bz2$/);
open P2, "$pair2_h" or die "cannot open $pair2\n";



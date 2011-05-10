#!/usr/bin/perl

=head1 NAME

pslx2sam.pl

=head1 DESCRIPTION

Convert a Blat output (modified PSLX format) into a SAM compatible output.

=head1 USAGE

pslx2sam.pl [OPTIONS]

OPTIONS
    Parameters         Description              Values      Default
    -i --input         Input file*              FILE        STDIN
    -o --output        Output file              FILE        STDOUT
    -u --unmap         Keep unmaped reads                   No
    -m --max-hits      Report up to Max hits    INT         10
    -p --poly-map      Report polymapped reads              No
    -h --help          Print this screen
    -v --verbose       Verbose mode
    --version          Print version
  
    * File can be compressed (gzip/bzip2)
  
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

# Main variables
my $version = 0.01;
my $input_h;
my ($sid, $seq, $nhit, $hits);
my ($mat, $dir, $cig, $chr, $pos);

# Parameters initialization
my $getversion = undef;
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $unmap      = undef;
my $polymap    = undef;
my $maxhits    = 10;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'u|unmap'         => \$unmap,
    'm|max-hits:i'    => \$maxhits,
    'p|polymap'       => \$polymap,
    'version'         => \$getversion
);
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $getversion);

if (defined $input) {
    $input_h = $input;
    $input_h = "gunzip  -c $input | " if ($input =~ m/\.gz$/);
    $input_h = "bunzip2 -c $input | " if ($input =~ m/\.bz2$/);
    open STDIN, "$input_h" or die "cannot open $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "cannot open $output\n";
}

while (<>) {
    chomp;
    ($sid, $seq, $nhit, $hits) = split (/\t/, $_);
    if ($nhit == 0) {
        printUnmap(\$sid, \$seq, \$nhits) if (defined $unmap);
    }
    if ($hits =~ m/_hap/) {
         filterHap(\$nhit, \$hits);
    }

    if ($nhits > $maxhits) {
        printUnmap(\$sid, \$seq, \$nhits) if (defined $polymap);
    }
    else {
        printSAM(\$sid, \$nhit, \$hits);
    }
    
}


#####################################################################
#
#                  S U B R O U T I N E S
#
#####################################################################

sub printVersion {
    print "pslx2sam.pl version $version\n";
    exit 1;
}
    

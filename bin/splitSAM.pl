#!/usr/bin/perl

=head1 NAME

splitSAM.pl

=head1 DESCRIPTION

Takes a SAM file/stream and wrote separate files for:
1) uniqelly mapped reads
2) unmapped reads
3) mapped reads in less than N hits (N != 1)
4) mapped reads in more than N hits

=head1 USAGE

perl splitSAM.pl -i SAM -1 UNIQ -u UNMAP -n N_HITS -m MAX_HITS -l 5

OPTIONS:
   Parameter         Description                    Values       Default
   -i  --input       Input SAM file* or stream      FILE         STDIN
   -1  --uniq        Uniquely mapped reads          FILE         uniqmap.sam
   -u  --unmap       No mapped reads                FILE         unmap.sam
   -n  --nmap        Reads with 1 < N < MAX hits    FILE         nmap.sam
   -p  --poly        Reads with > MAX hits          FILE         polymap.sam
   -l  --limit       Number of hits to define MAX   FILE         10
   -h  --help        Print help
   -v  --verbose     Verbose mode on
   --version         Print version

* File can be compressed (.gz|.bz2)
   
=head1 EXAMPLES

perl splitSAM.pl -i file.sam.gz

perl splitSAM.pl -l 5 < file.sam

zcat file.sam.gz | perl splitSAM.pl -1 my_uniq_map.sam

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

# Default parameters
my $our_version = 0.1;        # Version number
my $version  = undef;         # Version number call
my $help     = undef;         # Print help
my $verbose  = undef;         # Verbose mode
my $input    = undef;         # Input file or stream
my $uniq     = 'uniqmap.sam'; # Uniquely mapped reads output
my $unmap    = 'unmap.sam';   # Unmapped reads output
my $poly     = 'polymap.sam'; # Poly mapped reads output
my $nmap     = 'nmap.sam';    # N-mapped reads output
my $lim      = 10;            # Maximal hits to define $poly and $nmap

# Global variables
my $num_uniq  = 0;
my $num_poly  = 0;
my $num_unmap = 0;
my $num_nmap  = 0;
my $start     = time;
my $time      = undef;
my $end       = undef;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    '1|uniq:s'         => \$uniq,
    'u|unmap:s'        => \$unmap,
    'p|poly:s'         => \$poly,
    'n|nmap:s'         => \$nmap,
    'l|limit:i'        => \$lim,
    'version'          => \$version
    ) or pod2usage(-verbose => 2);
pod2usage(-verbose => 2) if (defined $help);

printVersion() if (defined $version);

# opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
open M, ">$uniq"  or die "cannot write file $uniq\n";
open U, ">$unmap" or die "cannot write file $unmap\n";
open N, ">$nmap"  or die "cannot write file $nmap\n";
open P, ">$poly"  or die "cannot write file $poly\n";

# Main operation
while (<>) {
    next if (m/^\@/); # skip headers
    my @line = split (/\t/, $_);
    if ($line[1] == 4) {
        $num_unmap++;
        print U $_;
    }
    else {
        if (m/NM:i:(\d+)/) {
            my $hits = $1;
            if ($hits == 0) {
                $num_uniq++;
                print M $_;
            }
            elsif ($hits < $lim) {
                $num_nmap++;
                print N $_;
            }
            else {
                $num_poly++;
                print P $_;
            }
        }
        else {
            die "NM: field is missing, aborting!\n";
        }
    }
}

$end = time;
$time = $end - $start;

if (defined $verbose) {
    warn "Maximal hits allowed  = $lim\n";
    warn "Unmapped reads        = $num_unmap\n";   
    warn "Uniquely mapped reads = $num_uniq\n";
    warn "N-mapped reads        = $num_nmap\n";
    warn "Polymapped reads      = $num_poly\n";
    warn "Time required (s)     = $time\n";
}

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}


#!/usr/bin/perl

=head1 NAME

countReads.pl

=head1 DESCRIPTION

Read a BAM file and count reads falling in a range set.

=head1 USAGE
  
  perl countReads.pl -b BAM -g GENES -o OUTPUT

  PARAMETER        DESCRIPTION                VALUE       DEFAULT
  -b --bam         BAM file                   FILE
  -g --genes       Gene annotation (BED)      FILE
  -o --out         Output file                FILE        STDOUT
  -u --uniq        Count only unique reads                No
  
  -h --help        Print this screen
  -v --verbose     Verbose mode
  

=head1 EXAMPLES

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
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help     = undef;      # Print help
my $verbose  = undef;      # Verbose mode
my $version  = undef;      # Version call flag
my $bam      = undef;      # BAM file
my $genes    = undef;      # Gene BED file
my $out      = undef;      # Output file
my $uniq     = undef;      # Count unique reads flag
# Define where is Samtools, adjust if required
#my $samtools = 'samtools view';   # samtools is in PATH
my $samtools = '/proj/hoodlab/share/programs/samtools/samtools view';

# Main variables
my $our_version = 0.1;     # Script version number
my %genes       = ();      # Gene coordinates
my %counts      = ();      # Reads per gene
my %seen        = ();      # Reads already counted (for spliced reads per chr)
my @reads       = ();      # Reads in a range
my ($chr, $ini, $end, $gen, $range, $reads, $read, $rid);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'g|genes:s'        => \$genes,
    'b|bam:s'          => \$bam,
    'u|uniq'           => \$uniq,
    'o|out:s'          => \$out
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $bam and defined $genes);
printVersion() if (defined $version);

if (defined $out) {
    open STDOUT, ">$out" or die "cannot write $out\n";
}

warn "loading gene coordinates from $genes\n" if (defined $verbose);
open G, "$genes" or die "cannot open $genes\n";
while (<G>) {
    chomp;
    ($chr, $ini, $end, $gen) = split (/\t/, $_);
    $ini++;
    push @{ $genes{$chr} }, "$ini,$end,$gen";
    $counts{$gen} = 0;
}
close G;

warn "querying the BAM file $bam\n" if (defined $verbose);
foreach $chr (keys %genes) {
    warn "    ... $chr\n";
    foreach $range (@{ $genes{$chr} }) {
        ($ini, $end, $gen) = split (/,/, $range);
        $reads = `$samtools $bam $chr:$ini-$end`;
        next unless (defined $reads);
        @reads = split (/\n/, $reads);
        foreach $read (@reads) {
            $read =~ m/NH:i:(\d+)/;
            $nh   = $1;
            if (defined $uniq) {
                next if ($nh > 1);
            }
            ($rid) = split (/\t/, $read);
            next if (defined $seen{$rid});
            $counts{$gen} += int(1 / $nh);
            $seen{$rid}    = 1;
        }
    }
    %seen = ();
}

warn "writing counts table\n" if (defined $verbose);
foreach $gen (sort keys %counts) {
    $reads = $counts{$gen};
    print "$gen\t$reads\n";
}

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}


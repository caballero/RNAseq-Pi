#!/usr/bin/perl

=head1 NAME

sam2coverage.pl

=head1 DESCRIPTION

Reads a SAM files (or BAM stream) and call the coverage per base.
The output is a text tab-delimited table with columns:
  1. chromosome/sequence
  2. position (1-based)
  3. coverage of the position

=head1 USAGE

sam2coverage.pl < SAM > COV

OPTIONS
    Parameter       Description                     Values        Default
    -s --sam        SAM file input*                 FILE          STDIN
    -o --out        Coverage output                 FILE          STDOUT
    -w --weight     Use weighted coverage**                       No
    -u --uniq       Use uniquelly mapped reads                    No
    -h --help       Print this screen
    -v --verbose    Verbose mode
    
     * File could be compressed (gzip|bunzip2)
    ** The weigthed mode use 1/N as a value for each position, where N is
       the total hits for the read, p.e. a read mapped in 2 different
       locations will support 0.5 in each base.

=head1 EXAMPLES

    # similar to a "pileup"
    sam2coverage.pl < map.sam > coverage.out 
    
    # weigthed mode
    sam2coverage.pl -w < map.sam > coverage.out 
    
    # all options
    sam2coverage.pl -w -s map.sam -o coverage.out -v
    
    # stream a BAM
    samtools view map.bam | sam2coverage.pl -w > coverage.out
    
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

# Global variables
my ($chr, $pos, $val, $cigar, $nhit, $sam_h, $out_h, $warn1);
my %data = ();

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $sam_f      = undef;
my $out_f      = undef;
my $weight     = undef;
my $uniq       = undef;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    's|sam:s'         => \$sam_f,
    'o|out:s'         => \$out_f,
    'w|weight'        => \$weight,
    'u|uniq'          => \$uniq
);
pod2usage(-verbose => 2) if (defined $help);


# Open files if required
if (defined $sam_f) {
    $sam_h = $sam_f;
    $sam_h = "gunzip  -c $sam_f | " if ($sam_f =~ m/gz$/);
    $sam_h = "bunzip2 -c $sam_f | " if ($sam_f =~ m/bz2$/);
    open STDIN, "$sam_h" or die "cannot open $sam_f\n";
}

if (defined $out_f) {
    open STDOUT, ">>$out_f" or die "cannot open $out_f\n";
}

# Parsing the SAM
warn "loading positions from SAM\n" if (defined $verbose);
while (<>) {
    next if (m/^\@/); # skip SAM headers
    if (m/NH:i:(\d+)/) {
        $nhit = $1;
    }
    else {
        $warn1++;
        warn "NH variable not found, skiping weights\n" if (defined $verbose and $warn1 == 1);
        $nhit = 1;
    }
    
    if (defined $uniq) {
        next unless ($nhit == 1);
    }
    
    $nhit = 1 unless (defined $weight);
    my @line = split (/\t/, $_);
    next if ($line[1] == 4); # skip unmapped reads
    $chr   = $line[2];
    $pos   = $line[3];
    $cigar = $line[5];
    $val   = 1 / $nhit;
    # decode CIGAR
    while ($cigar =~ m/(\d+)(\w)/g) {
        my $num = $1;
        my $cig = $2;
        if ($cig =~ m/[SIDH]/) {
            next;
        }
        elsif ($cig eq 'N') {
            $data{$chr}{$pos - 1}{'s'}++;
            $pos += $num;
            $data{$chr}{$pos}{'s'}++;
        }
        elsif ($cig eq 'M') {
            for (my $i = 1; $i <= $num; $i++) {
                $data{$chr}{$pos}{'w'} += $val;
                if ($nhit == 1) {
                    $data{$chr}{$pos}{'u'}++;
                }
                else {
                    $data{$chr}{$pos}{'p'}++;
                }
                $pos++;
            }
        }
    }
}

# Printing coverages per base
warn "printing coverages\n" if (defined $verbose);
foreach $chr (keys %data) {
    foreach $pos (sort {$a<=>$b} keys %{ $data{$chr} }) {
        my $w = 0; $w = sprintf("%.2f", $data{$chr}{$pos}{'w'}) if (defined $data{$chr}{$pos}{'w'});
        my $u = 0; $u = $data{$chr}{$pos}{'u'} if (defined $data{$chr}{$pos}{'u'});
        my $p = 0; $p = $data{$chr}{$pos}{'p'} if (defined $data{$chr}{$pos}{'p'});
        my $s = 0; $s = $data{$chr}{$pos}{'s'} if (defined $data{$chr}{$pos}{'s'});
        print "$chr\t$pos\t$w\t$u\t$p\t$s\n";
    }
}


warn "done\n" if (defined $verbose);

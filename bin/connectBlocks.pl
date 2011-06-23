#!/usr/bin/perl

=head1 NAME

connectBlocks.pl

=head1 DESCRIPTION

Report the number of reads "connecting" (splices) the detected "blocks" (exons).
The inputs are the output of cov2block.pl and the SAM file.
Output is a table showing (text tab-delimited):
1. Connection ID
2. Block 1
3. Block 2
4. Reads

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -s --sam         SAM file                   File*      STDIN
    -i --input       Blocks input               File*       
    -o --output      Output                     File       STDOUT
    -m --min         Minimal reads supporting   INT        1
    -l --label       Laber connection           STR        CONN
    -h --help        Print this screen
    -v --verbose     Verbose mode
    
    * File can be compressed (gzip|bzip2)

=head1 EXAMPLES

   connectBlocks.pl -s sample.sam -i sample.block -o sample.conn

   connectBlocks.pl -s sample.sam -i sample.block -o sample.conn -m 3
   
   connectBlocks.pl -s sample.sam.gz -i sample.block -o sample.conn -l SAMPLE

   samtools view sample.bam | connectBlocks.pl -i sample.block > sample.conn

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
my ($id, $chr, $pos, $cigar, $ini, $end);
my ($blk, $cnn);      # counters
my ($sam_h $input_h); # filehandlers
my %blocks = ();

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $sam        = undef;
my $min        = 1;
my $label      = 'CONN';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    's|sam:s'         => \$sam,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'm|min:i'         => \$min,
    'l|label:s'       => \$label
);
pod2usage(-verbose => 2)     if (defined $help);
pod2usage(-verbose => 2) unless (defined $input);


# Opening files
$input_h = $input;
$input_h = "gunzip  -c $input | " if ($input =~ m/gz$/);
$input_h = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
open IN, "$input_h" or die "cannot open $input\n";

if (defined $sam) {
    $sam_h = $sam;
    $sam_h = "gunzip  -c $sam | " if ($sam =~ m/gz$/);
    $sam_h = "bunzip2 -c $sam | " if ($sam =~ m/bz2$/);
    open STDIN, "$sam_h" or die "Cannot read file $sam\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

# reading blocks
warn "loading blocks\n" if (defined $verbose);
while (<IN>) {
    chomp;
    ($id, $chr, $ini, $end) = split (/\t/, $_);
    $blocks{$chr}{$ini}{$end} = $id;
    $blk++;
}
warn "$blk blocks found\n" if (defined $verbose);

while (<>) {
    my @line = split (/\t/, $_);
    $chr   = $line[2];
    $pos   = $line[3];
    $cigar = $line[5];
    next unless ($cigar =~ m/N/); # is a spliced read?
    # define borders
    my @tag = split (/(\d+\w)/, $_);
    foreach my $tag (@tag) {
        if ($tag =~ m/M$/) {
            chop $tag;
            $pos += $tag;
        }
        elsif ($tag =~ m/N$/) {
            chop $tag;
            $left   = $pos;
            $pos   += $tag;
            $rigth  = $pos;
            $splice = checkBlocks($chr, $left, $right);
            $cnn++ if ($splice == 1);
        }
        else {
            # do nothing
        }
    }
}

sub checkBlocks {
    my ($c, $l, $r) = @_;
    my $res = 0;
    my $b1  = undef;
    my $b2  = undef;
    foreach my $i (keys %{ $blocks{$c} }) {
        
    }
    return $res;
}

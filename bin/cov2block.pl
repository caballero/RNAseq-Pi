#!/usr/bin/perl

=head1 NAME

cov2block.pl

=head1 DESCRIPTION

Reads a coverage file/stream and detect expression blocks (exons). 

The coverages file is a text tab-delimited file with columns:
1. sequence/chromosome id
2. position (1-based)
3. coverage per base
This file must be sorted by columns 1-2 (sort -k1,1 -kn2,2)

The output is a text tab-delimited file with columns:
1. sequence/chromosome id
2. start position
3. end position
4. average coverage

=head1 USAGE

cov2block.pl < COVERAGES > BLOCKS

OPIONS
    Parameter        Description                  Value       Default
    -i --input       Coverages*                   File        STDIN
    -o --output      Output blocks                File        STDOUT
    -w --window      Slicing window size          INT         100
    -s --step        Overlap between windows      INT         20
    -m --mincov      Minimal coverage             FLOAT       1.0
    
    -h --help        Print this screen
    -v --verbose     Verbose mode
    
    * File can be compressed (gzip|bzip2)

=head1 EXAMPLES

   cov2block.pl < maped.cov > detected.block

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

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $window     =   100;
my $step       =    20;
my $mincov     =   1.0;

# Global variables
my ($seq, $pos, $cov, $input_h);
my %block = ();

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|ouput:s'       => \$output,
    'w|window:i'      => \$window,
    's|step:i'        => \$step,
    'm|mincov:s'      => \$mincov
);
pod2usage(-verbose => 2) if (defined $help);

# Opening files if required
if (defined $input) {
    warn "Reading input from: $input\n" if (defined $verbose);
    $input_h = $input;
    $input_h = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input_h = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input_h" or die "cannot open $input\n";
}

if (defined $output) {
    warn "Writing output in:  $output\n" if (defined $verbose);
    open STDOUT, ">$output" or die "cannot open $output\n";
}
warn "Window size       = $window\n" if (defined $verbose);
warn "Step size         = $step\n" if (defined $verbose);
warn "Minimal coverage  = $mincov\n" if (defined $verbose);

while (<>) {
    chomp;
    my ($seq, $pos, $cov) = split (/\t/, $_);
    next if ($cov < $mincov);
    if (defined $block{'seq'}) { #block is being expanded
        if ($seq ne $block{'seq'}) {
            warn "processed sequence: $block{'seq'}\n" if (defined $verbose);
            printBlock();
            newBlock($seq, $pos, $cov);
        }
        elsif ($pos > $block{'end'} + $window) {
            printBlock();
            newBlock($seq, $pos, $cov);
        }
        elsif ($block{'num'} >= $step) {
            $block{'mean'} = mean(@{ $block{'cov'} });
            $block{'sd'}   =   sd($block{'mean'}, @{ $block{'cov'} });
           
            if( $cov >= ($block{'mean'} - (2*$block{'sd'}) - 1) and $cov <= ($block{'mean'} + (2*$block{'sd'}) + 1)) {
                $block{'end'} = $pos;
                $block{'num'}++;
                push @{ $block{'cov'} }, $cov;
            }
            else {
                printBlock();
                newBlock($seq, $pos, $cov);
            }
        }
        else {
            $block{'end'} = $pos;
            $block{'num'}++;
            push @{ $block{'cov'} }, $cov;
        }   
    }
    else { # a new block
        newBlock($seq, $pos, $cov);
    }
}
printBlock() if (defined $block{'num'}); # last block call
warn "processed sequence: $block{'seq'}\n" if (defined $verbose);


################################################################
#           S  U  B  R  O  U  T  I  N  E  S                    #
################################################################

sub newBlock {
    my ($s, $p, $c) = @_;
    %block = ();
    $block{'seq'}    = $s;
    $block{'ini'}    = $p;
    $block{'end'}    = $p;
    $block{'cov'}[0] = $c;
    $block{'num'}    = 1;
}

sub printBlock {
    if($block{'num'} >= $step) {
        $block{'mean'} = sprintf("%.2f", mean(@{ $block{'cov'} }));
        $block{'sd'}   = sprintf("%.2f", sd($block{'mean'}, @{ $block{'cov'} }));
        print join ("\t", $block{'seq'}, $block{'ini'}, $block{'end'}, $block{'mean'}, $block{'sd'});
        print "\n";
    }
}

sub mean {
    my $sum  = 0;
    my $num  = 0;
    foreach my $x (@_) {
        $sum += $x;
        $num++;
    }
    if ($num > 0) {
        return $sum / $num;
    }
    else {
        return 'NA';
    }
}

sub sd {
    my $mean = shift @_;
    my $sum  = 0;
    my $num  = 0;
    foreach my $x (@_) {
        $sum += ($x - $mean) * ($x - $mean);
        $num++;
    }
    if ($num > 1) {
        return sqrt($sum / ($num - 1));
    }
    else {
        return 'NA';
    }
}

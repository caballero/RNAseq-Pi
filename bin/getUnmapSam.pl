#!/usr/bin/perl

=head1 NAME

parseSAM.pl

=head1 DESCRIPTION

Read a SAM file input and separate mapped and unmapped reads (mapped reads are
reported in SAM format, unmapped reads can be reported as Fasta/Fastq).

=head1 USAGE

perl parseSAM.pl [OPTIONS] < SAM > OUTPUT

OPTIONS:
   Parameter        Description               Values          Default
   -i --input       Input SAM file*           File name       STDIN
   -o --output      Output file               File name       STDOUT
   -f --format      Output format             fa/fq           fq
   -m --mapped      Write mapped reads here   File name
   -q --qual        Use this quality score    QS              B
   -v --verbose     Verbose mode
   -h --help        Print this screen
   
   * File can be compressed (gzip/bzip2) but requires option -i
   
=head1 EXAMPLES

perl parseSAM.pl < align.sam > unmapped_reads.fq

perl parseSAM.pl -i align.sam.gz -f fa -m mapped.sam -o unmapped_reads.fa

samtools view align.bam | perl parseSAM.pl > unmmaped_reads.fq

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

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help    = undef;    # print help
my $verbose = undef;    # verbose mode
my $format  = 'fq';     # sequence format
my $input   = undef;    # input file
my $output  = undef;    # output file
my $mapped  = undef;    # report mapped reads here
my $dqual   = 'B';      # default quality score

# Main variables
my $num_mapped   = 0;
my $num_unmapped = 0;
my $total        = 0;

# Calling options
GetOptions(
    'h|help'      => \$help,
    'v|verbose'   => \$verbose,
    'i|input:s'   => \$input,
    'f|format:s'  => \$format,
    'o|output:s'  => \$output,
    'm|mapped:s'  => \$mapped,
    'q|qual:s'    => \$dqual
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if (defined $help);

# Opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}
if (defined $mapped) {
    open MAPPED, ">$mapped" or die "cannot write file $mapped\n";
}

# Reading SAM file
while (<>) {
    # catch SAM header
    if (/^\@/) {
        #print MAPPED $_ if (defined $mapped);
    }
    else {
        my @line = split (/\t/, $_);
        my $flag = $line[1];
        
        if ($flag == 4) { # unmapped read
            my $id   = $line[0];
            my $seq  = $line[9];
            my $qual = $line[10];
            printSeq($id, $seq, $qual);
            $num_unmapped++;
        }
        else { # mapped read
            print MAPPED $_ if (defined $mapped);
            $num_mapped++;
        }
    }
}

if (defined $verbose) {
    my $total     = $num_mapped + $num_unmapped;
    my $per_map   = sprintf( "%.2f", 100 * $num_mapped   / $total);
    my $per_unmap = sprintf( "%.2f", 100 * $num_unmapped / $total);
    warn <<__INFO__    
# Total sequences    = $total
# Mapped sequences   = $num_mapped ($per_map\%)
# Unmapped sequences = $num_unmapped ($per_unmap\%)
__INFO__
;
}

###################################
####   S U B R O U T I N E S   ####
###################################

# printSeq -> print a Fasta/Fastq sequence
# call: printSeq( $seq_id, $seq_nuc, $seq_qual )
# return: nothing

sub printSeq {
    my ($i, $s, $q) = @_;
    if ($format eq 'fa') {
        print ">$i\n$s\n";
    }
    elsif ($format eq 'fq') {
        if ($q eq '*') {
            $q = $dqual x length $s;
        }
        print "\@$i\n$s\n+\n$q\n";
    }
    else {
        die "format is not recognized: $format\n";
    }
}

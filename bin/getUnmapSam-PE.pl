#!/usr/bin/perl

=head1 NAME

getUnmapSam-PE.pl

=head1 DESCRIPTION

Read a SAM file input and separate mapped and unmapped reads (mapped reads are
reported in SAM format, unmapped reads can be reported as Fasta/Fastq).

Special version for paired-end sequences.

=head1 USAGE

perl parseSAM.pl [OPTIONS] < SAM > OUTPUT

OPTIONS:
   Parameter        Description               Values          Default
   -i --input       Input SAM file*           File name       STDIN
   -p --prefix      Output file prefix**      File name       
   -1 --out1        Output file pair 1        File name
   -2 --out2        Output file pair 2        File name
   -f --format      Output format             fa/fq           fq
   -m --mapped      Write mapped reads here   File name
   -q --qual        Use this quality score    QS              B
   -v --verbose     Verbose mode
   -h --help        Print this screen
   
   *  File can be compressed (gzip/bzip2) but requires option -i
   ** Prefix value will be used to create 2 output files: PRE-1.fq, PRE-2.fq
   
=head1 EXAMPLES

perl parseSAM.pl -p unmap_pair < align.sam

perl parseSAM.pl -i align.sam.gz -f fa -m mapped.sam -1 unmap_1.fa -2 unmap_2.fa

samtools view align.bam | perl parseSAM.pl -p unmapped_reads

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
my $prefix  = undef;    # prefix name for outputs
my $out1    = undef;    # output 1
my $out2    = undef;    # output 2
my $mapped  = undef;    # report mapped reads here
my $dqual   = 'B';      # default quality score

# Main variables
my $num_mapped   = 0;
my $num_unmap1   = 0;
my $num_unmap2   = 0;

# Calling options
GetOptions(
    'h|help'      => \$help,
    'v|verbose'   => \$verbose,
    'i|input:s'   => \$input,
    'f|format:s'  => \$format,
    'p|prefix:s'  => \$prefix,
    '1|out1:s'    => \$out1,
    '2|out2:s'    => \$out2,
    'm|mapped:s'  => \$mapped,
    'q|qual:s'    => \$dqual
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $prefix or (defined $out1 and defined $out2));

# Opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $mapped) {
    open MAPPED, ">$mapped" or die "cannot write file $mapped\n";
}

if (defined $prefix) {
    $out1 = "$prefix-1.$format";
    $out2 = "$prefix-2.$format";
}
open O1, ">$out1" or die "cannot write file $out1\n";
open O2, ">$out2" or die "cannot write file $out2\n";


# Reading SAM file
while (<>) {
    # catch SAM header
    if (/^\@/) {
        #print MAPPED $_ if (defined $mapped);
    }
    else {
        my @line = split (/\t/, $_);
        my $flag = $line[1];
        my $id   = $line[0];
        my $seq  = $line[9];
        my $qual = $line[10];
        my $fseq = undef;
        
        if ($flag == 77) { # unmapped read 1
            $id =~ '/1';
            $fseq = formatSeq($id, $seq, $qual);
            print O1 $fseq;
            $num_unmap1++;
        }
        elsif ($flag == 141) { # unmapped read 2
            $id .=~ '/2';
            $fseq = formatSeq($id, $seq, $qual);
            print O2 $fseq;
            $num_unmap2++;
        }
        else { # mapped read
            print MAPPED $_ if (defined $mapped);
            $num_mapped++;
        }
    }
}

if (defined $verbose) {
    warn <<__INFO__    
# Mapped sequences = $num_mapped
# Unmapped seqs #1 = $num_unmap1
# Unmapped seqs #2 = $num_unmap1
__INFO__
;
}

###################################
####   S U B R O U T I N E S   ####
###################################

# printSeq -> print a Fasta/Fastq sequence
# call: printSeq( $seq_id, $seq_nuc, $seq_qual )
# return: nothing

sub formatSeq {
    my ($i, $s, $q) = @_;
    my $r = undef;
    if ($format eq 'fa') {
        $r = ">$i\n$s\n";
    }
    elsif ($format eq 'fq') {
        if ($q eq '*') {
            $q = $dqual x length $s;
        }
        $r = "\@$i\n$s\n+\n$q\n";
    }
    else {
        die "format is not recognized: $format\n";
    }
    return $r;
}

#!/usr/bin/perl -w
use strict;

=head1 NAME

getUnampSam.pl 

=head1 DESCRIPTION

Read a SAM file, get the unmapped reads, write a Fastq/Fasta file with them.

=head1 USAGE

perl getUnmapSam.pl < SAM > FASTQ

OPTIONS:
   Parameter       Description                  Values     Default
   -i --input      Read sequences from here     File*      STDIN
   -o --output     Write sequences here         File       STDOUT
   -f --format     Output format (Fastq/Fasta)  fq/fa      fq
   -g --group      Write sequences in batches   int        100000
   -v --verbose    Verbose mode
   -h --help       Print help
   
   *File can be compressed (gzip/bzip2) but requires option -i
   
=head1 EXAMPLES

perl getUnmapSam.pl < SAM > FASTQ

perl getUnmapSam.pl -f fa -i SAM -o FASTA

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
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help     = undef;         # print help
my $verbose  = undef;         # verbose mode
my $format   = 'fq';          # sequence format
my $input    = undef;         # input SAM file
my $output   = undef;         # output file
my $batch    = 100000;        # size of the batch

# Main variables
my $id       = undef;         # Sequence ID
my $bit      = undef;         # SAM bit flag
my $seq      = undef;         # Sequence
my $qual     = undef;         # Fastq qualities
my %seqs     = ();            # Hash to keep sequences in memory
$seqs{'num'} = 0;             # Sequences counter
my $nbatch   = 0;             # Counter for batches
my @arr      = ();            # array for data

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output,
    'f|format:s'       => \$format,
    'b|batch:i'        => \$batch,
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);

# opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}

while (<>) {
    next if (m/^\@/); # skip SAM headers
    @arr  = split (/\t/, $_);
    $id   = $arr[0];
    next if (defined $seqs{'used'}{$id}); # No duplicates
    $bit  = $arr[1];
    $seq  = $arr[8];
    $qual = $arr[9];
    next unless ($bit == 4); # Only unmapped reads (flag = 4)
    $seqs{'num'}++;
    $seqs{'used'}{$id} = 1;
    
    if ($format eq 'fq') {
        $seqs{'seqs'} .= "\@$id\n$seq\n+\n$qual\n";
    }
    elsif ($format eq 'fa') {
        $seqs{'seqs'} .= ">$id\n$seq\n";
    }
    else {
        die "Sorry, I don't recognize this format: $format\n";
    }
    
    flushSeqs() if ($seqs{'num'} >= $batch);
}
flushSeqs();

# Print report in verbose mode 
warn "found $nbatch sequences\n" if (defined $verbose);


###################################
####   S U B R O U T I N E S   ####
###################################

# flushSeqs -> print sequences
# Call: flushSeqs()
# Return: nothing
sub flushSeqs  {
    $nbatch += $seqs{'num'};
    warn "processed $nbatch seqs\n" if (defined $verbose);
    print $seqs{'seqs'} if (defined $seqs{'seqs'});
    $seqs{'num'}  = 0;
    $seqs{'seqs'} = '';
}


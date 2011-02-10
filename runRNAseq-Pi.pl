#!/usr/bin/perl

=head1 NAME

runRNAseq-Pi.pl

=head1 DESCRIPTION

Main script to run RNAseq-Pi.

=head1 USAGE

perl runRNAseq-Pi.pl -c config.txt -i READS -o SAM.OUT [OPTIONS]

OPTIONS:
   Parameter       Description            Values         Default
   -i --input      Read files             File names*    None, required
   -o --output     Final SAM output       File name      None, required
   -c --config     Configuration file     File name**    None, required
   -f --format     Sequences format       fq/fa          fq
   -t --threads    Multicore number       Int***         1
   -b --batch      Analyze seqs in batch  Int            1000000
   -p --paired     Activate Paired-End                   No
   -v --verbose    Verbose mode                          No
   -h --help       Print this screen                     No
   
   Notes:
   *   Files can be fasta/fastq, compressed (gzip or bzip2), multiple files can
       be separated with ",", but paired-sequences must be grouped with ":".
   **  See "config.txt" for a complete description.
   *** Launch job using N processors (bwa and bowtie can use them).
   
=head1 EXAMPLES

perl runRNAseq-Pi.pl

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
use RNAseqPi;
use Getopt::Long;
use Pod::Usage;

# Global variables
my %param    = ();
my $input    = undef;
my $output   = undef;
my $config   = undef;
my $batch    = 1000000;
my $help     = undef;
my $verbose  = undef;
my $threads  = 1;
my $paired   = undef;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input=s'        => \$input,
    'o|output=s'       => \$output,
    'c|config=s'       => \$config,
    'b|batch:i'        => \$batch,
    't|thread:i'       => \$threads,
    'p|paired'         => \$paired
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $input and defined $config and defined $output);


loadParam();



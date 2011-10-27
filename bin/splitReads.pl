#!/usr/bin/perl

=head1 NAME

splitReads.pl

=head1 DESCRIPTION

Take a fasta/fastq file/stream and divide it into N files.

=head1 USAGE

OPTIONS
    Parameter        Description                 Value      Default
    -i --input       Input*                      File       STDIN
    -p --prefix      Output prefix               NAME       reads
    -n --num         Split the file in N         INT        4
    -f --format      Input format (fastq/fasta)  fq/fa      fq
    
    -h --help        Print this screen
    -v --verbose     Verbose mode
    
    * Input can be compressed (gzip/bzip2) but need -i

=head1 EXAMPLES

    - Split a fasta file in 2 parts:
    splitReads.pl -f fa -n 2 -i file.fa.gz -p file_part.

    - Split a fastq file in 8 parts:
    splitReads.pl -n 8 -p file_part. < file.fq
   
    
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

#use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Global
my ($out, $line, $seq, $nl, $out_fh_max);
my $out_fh = 'fh00';

# Parameters initialization
my $help       =   undef;
my $verbose    =   undef;
my $input      =   undef;
my $prefix     = 'reads';
my $slices     =       4;
my $format     =    'fq';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'p|prefix:s'      => \$prefix,
    'n|num:i'         => \$slices,
    'f|format:s'      => \$format
);
pod2usage(-verbose => 2) if (defined $help);

# define how to read the sequences
if    ($format eq 'fq') { 
    $line = 4; 
}
elsif ($format eq 'fa') { 
    $line = 2; 
}
else { 
    die "Format not recognized: $format\n"; 
}

# open input if required
if (defined $input) {
    my $input_h = $input;
    $input_h = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input_h = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input_h" or die "Cannot read file $input\n";
}

# create outputs
$out = $out_fh;
for (my $i = 1; $i <= $slices; $i++) { 
    $out++;
    open $out, ">$prefix.$i" or die "cannot open $prefix.$i\n";
}
$out_fh_max = $out;

# read input
$nl  = 0;
$out = $out_fh;
while (<>) {
    $nl++;
    $seq .= $_;
    if ($nl == $line) {
        $out++;
        print $out $seq;
        $out = $out_fh if ($out eq $out_fh_max);
        $nl  = 0;
        $seq = '';
    }
}

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

sampleReads.pl

=head1 USAGE

sampleReads.pl -n N -t TOTAL -i FASTQ -o OUTPUT -f fq -g fa

OPTIONS:
    Parameter         Description                       Value       Default
    -n  --num         Total reads to select             Int
    -t  --total       Total reads in the dataset        Int
    -i  --in          Input file                        File
    -o  --out         Output file                       File
    -f  --format_in   Format of infile (fasta/fastq)    fa/fq       fq
    -g  --format_out  Format of outfile (fasta/fastq)   fa/fq       fq
    -q  --quality     Use this quality score (fa->fq)   [2-40]      9
    -h  --help        Print this screen

=head1 DESCRIPTION

Randomly choose a sample of N reads.

=cut

# Variables 
my $help       = undef;
my $num        = undef;
my $tot        = undef;
my $input      = undef;
my $output     = undef;
my $format_in  = 'fq';
my $format_out = 'fq';
my $qual       = 9;
my %select     = ();
my $comm       = undef;
my $line       = 0;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output,
    'f|format_in:s'    => \$format_in,
    'g|format_out:s'   => \$format_out,
    'q|quality:i'      => \$qual,
    'n|num=i'          => \$num,
    't|tot=i'          => \$tot
) or pod2usage(-verbose => 2);
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $num and defined $tot);

my $codeq = encodeQual($qual);
selectPositions();

# Opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}

if    ($format_in eq 'fq') { $comm = '@'; }
elsif ($format_in eq 'fa') { $comm = '>'; }
else  { die " cannot recognize this format $format_in\n"; }

$/ = "\n$comm";
while (<>) {
    $line++;
    next unless (defined $select{$line});
    s/$comm//;
    if    ($format_in eq $format_out)                  { print $comm . $_; }
    elsif ($format_in eq 'fq' and $format_out eq 'fa') { print fq2fa($_);  }
    elsif ($format_in eq 'fa' and $format_out eq 'fq') { print fa2fq($_);  }
    else {       die "not recognized, in: $format_in, out: $format_out\n"; }
}

# SUBROUTINES

sub selectPositions {
    for (my $i = 0; $i <= $num; $i++) {
        my $pos = int(rand $tot);
        if (defined $select{$pos}) {
            $i--;
        } else {
            $select{$pos} = 1;
        }
    }
}

sub fq2fa {
    my $fq = shift @_;
    my ($id, $seq, $sep, $qual) = split (/\n/, $fq);
    return ">$id\n$seq\n";
}

sub fa2fq {
    my $fa = shift @_;
    my ($id, $seq) = split (/\n/, $fa);
    my $qual = $codeq x length $seq;
    return "\@$id\n$seq\n+\n$qual\n"
}

sub encodeQual {
    my $q = shift @_;
    my $s = chr($q + 64);
    return $s;
}

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

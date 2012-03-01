#!/usr/bin/perl

=head1 NAME

correctBamQual.pl

=head1 DESCRIPTION

Change the quality scores in a SAM/BAM file using the values from a FASTQ file.

*warning this script load all reads in the Fastq to memory, so be sure you 
have enough.

=head1 USAGE

perl correctBamQual.pl -b BAM -f FASTQ -o OUTPUT

    PARAMETER           DESCRIPTION                 VALUE        DEFAULT
    -b --bam            BAM file                    FILE         
    -s --sam            SAM file                    FILE         STDIN  
    -o --out            Output (SAM file)           FILE         STDOUT
    -f --fastq          Fastq file                  FILE
    -c --convert        Convert SAM output to BAM
    
    -h --help           Print this screen
    -v --verbose        Verbose mode

=head1 EXAMPLES

    perl correctBamQual.pl -b ALIGN.bam -f READS.fq -o OUTPUT.sam
    
    perl correctBamQual.pl -s ALIGN.sam.gz -f READS.fq.gz -c -o OUTPUT.bam
    
    samtools view ALIGN.bam | perl correctBamQual.pl -f READS.fq > OUTPUT.sam

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
my $help     = undef;         # Print help
my $verbose  = undef;         # Verbose mode
my $version  = undef;         # Version call flag
my $bam      = undef;
my $sam      = undef;
my $fastq    = undef;
my $out      = undef;
my $convert  = undef;
#my $samtools = 'samtools';
my $samtools = '/proj/hoodlab/share/programs/samtools/samtools';

# Main variables
my $our_version = 0.1;        # Script version number
my %quals = ();
my ($read, $qual, $line);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose
) or pod2usage(-verbose => 2);
printVersion() if (defined $version);
   
pod2usage(-verbose => 2) if (defined $help);
die "please specify only one option for -s|-b\n" if (defined $bam and defined $sam);
die "please specify a FASTQ file (-f)\n"     unless (defined $fastq);
die "please specify a SAM output (-o)\n"         if (defined $convert and !defined $out);

# Open files (if required)
if (defined $bam) {
    my $bam_fh = defineFH($bam);
    open STDIN, "$bam_fh" or die "cannot open file $bam\n";
}
elsif (defined $sam) {
    my $sam_fh = defineFH($sam);
    open STDIN, "$sam_fh" or die "cannot open file $sam\n";
}

if (defined $out) {
    open STDOUT, ">$out" or die "cannot open file $out\n";
}

my $fastq_fh = defineFH($fastq);
open FQ, "$fastq_fh" or die "cannot open file $fastq_fh\n";

warn "loading quality scores\n" if (defined $verbose);
$line = 0;
while (<FQ>) {
    chomp;
    $line++;
    if ($line == 1) {
        s/^\@//;
        s/#\d+\/\d+$//;
        $read = $_;
    }
    elsif ($line == 4) {
        $quals{$read} = $_;
        $line = 0;
    }
}

warn "changing quality scores\n" if (defined $verbose);
while (<>) {
    my @line = split (/\t/, $_);
    $line[0] =~ s/#\d+\/\d+$//;
    if (defined $quals{ $line[0] }) {
        $line[10] = $quals{ $line[0] }; 
    }
    else {
        warn " .. missing quality scores for $line[0]\n" if (defined $verbose);
    }
    print join "\t", @line;
}

if (defined $convert) {
    warn "converting SAM to BAM\n" if (defined $verbose);
    my $bam_out = $out;
    $bam_out =~ s/sam$/bam/;
    system ("samtools view -bS $out > $bam_out");
    unlink $out;
}

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub defineFH {
    my $fi = shift @_;
    my $fh = $fi;
    if    ($fi =~ m/.gz$/) { $fh = "gunzip      -c $fi | "; }
    elsif ($fi =~ m/bz2$/) { $fh = "bunzip2     -c $fi | "; }
    elsif ($fi =~ m/bam$/) { $fh = "$samtools view $fi | "; }
    
    return $fh;
}

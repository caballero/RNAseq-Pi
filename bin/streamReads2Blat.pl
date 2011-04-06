#!/usr/bin/perl

=head1 NAME

streamReads2Blat.pl

=head1 DESCRIPTION


=head1 USAGE


=head1 EXAMPLES


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

# Default parameters
my $help      = undef;         # print help
my $verbose   = undef;         # verbose mode
my $version   = undef;         # version call flag
my $input     = undef;         # input file
my $format    = 'fq';          # input format
my $output    = undef;         # output file


# Main variables
my $index_dir = '/proj/hoodlab/share/programs/blat-indexes'; # Blat indexes dir
my $fasta     = 'read.fa';
my $psl       = 'read.psl';
my $gfclient  = '/proj/hoodlab/share/programs/blat/gfClient';
my $host    = 'localhost';
my @indexes   = (   'solexa_primers.2bit',
                    'human.GRCh37.61.rRNA-MT.2bit',
                    'human_RepBase15.10.2bit',
                    'ERCC_reference_081215.2bit',
                    'hs37.61.2bit'                     );
my %indexes  = ();
$indexes{'solexa_primers.2bit'}{'name'}          = 'primer match';
$indexes{'solexa_primers.2bit'}{'port'}          = 1111;
$indexes{'human.GRCh37.61.rRNA-MT.2bit'}{'name'} = 'rRNA/MT match';
$indexes{'human.GRCh37.61.rRNA-MT.2bit'}{'port'} = 1112;
$indexes{'human_RepBase15.10.2bit'}{'name'}      = 'repeat match';
$indexes{'human_RepBase15.10.2bit'}{'port'}      = 1113;
$indexes{'ERCC_reference_081215.2bit'}{'name'}   = 'ERCC macth';
$indexes{'ERCC_reference_081215.2bit'}{'port'}   = 1114;
$indexes{'hs37.61.2bit'}{'name'}                 = 'hg19 match';
$indexes{'hs37.61.2bit'}{'port'}                 = 1115;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'f|format:s'       => \$format,
    'o|output:s'       => \$output,
    'version'          => \$version
    ) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if(defined $version);

# Opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
    open STDOUT, ">$output" or die "cannot open $output\n";
}

# Main program
if ($format eq 'fq') {
    $/ = "\n\@";
    while (<>) {
        my ($id, $seq, $sep, $qual) = split (/\n/, $_);
        $id =~ s/\@//;
        writeFa($id, $seq);
        my $hit = searchHit($id, $seq);
        print $hit;
    }
    unlink "$fasta";
    unlink "$psl";
}
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
        my ($id, $seq) = split (/\n/, $_);
        $id =~ s/>//;
        writeFa($id, $seq);
        my $hit = searchHit($id, $seq);
        print $hit; 
    }
}
else {
    die "Format $format is not FASTQ/FASTA\n";
}

unlink "$fasta";
unlink "$psl";


# SUBROUTINES

sub writeFa {
    my ($id, $seq) = @_;
    open  FA, ">$fasta" or die "cannot open $fasta\n";
    print FA ">$id\n$seq\n";
    close FA;
}

sub searchHit {
    my ($id, $seq) = @_;
    my $res = "$id\t$seq\tNo_hit_found";
    foreach my $target (@indexes) {
        my $name = $indexes{$target}{'name'};
        runBlat($target);
        my $hit = checkHit();
        if (defined $hit) {
            $res = "$id\t$seq\t$name\t$hit";
            last;
        }
    }
    return $res;
}

sub runBlat {
    my $target = shift @_;
    my $port   = $indexes{$target}{'port'};
    system ("$gfclient $host $port . -nohead $fasta $psl");
}

sub checkHit {
    my $res = undef;
    if (-s $psl) {
        open PSL, "$psl" or die "cannot open $psl\n";
        my $best = -1;
        while (<PSL>) {
            my @array = split (/\t/, $_);
            my $score = $array[0] - $array[1];
            if ($score > $best) {
                $res = join ":", @array;
            }
            elsif ($score == $best) {
                $res .= '|';
                $res = join ":", @array;
            }
        }
        close PSL;
    }
    return $res;
}
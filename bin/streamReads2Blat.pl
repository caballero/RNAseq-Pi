#!/usr/bin/perl

=head1 NAME

streamReads2Blat.pl

=head1 DESCRIPTION

Create a pseudo stream to query a file with reads to a serial blat searches.
Each read is classified in one of the index o declared not found.
The Blat gfServers must be active with the same index/port configuration (use
blatServers.pl script before).

=head1 USAGE

perl streamReads2Blat.pl -i FASTQ -o ALIGN

OPTIONS
    Parameters         Description                         Value     Default
    -i --input         Input file*                         FILE      STDIN
	-o --output        Output file                         FILE      STDOUT
	-f --format        Input format                        fq/fa     fq
	-s --minscore      Minimal Blat score**                INT       100
	-p --percent       Minimal identity percent            INT       90
	-m --maxintron     Maximal intron size                 INT       750000
	-n --nreads        Report time every N reads           INT       1000
	-h --help          Print this screen
	-v --verbose       Activate verbose mode
	--version          Print version number
	
 * Input file can be compressed (.gz|.bz2)
** Blat score = (2 * NumMatches) - NumMismatches

=head1 EXAMPLES

perl streamReads2Blat.pl -i reads.fq.gz -o align.out

perl streamReads2Blat.pl -i reads.fa -f fa -o align.out -s 50 -p 95 -n 10000

perl streamReads2Blat.pl < reads.fq > align.out

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
my $minscore  = 100;           # minimal blat score
my $minident  = 90;            # minimal identity percent
my $block     = 1000;          # partial report every N reads;
my $maxintron = 750000;        # max Intron size

# Main variables
my $our_version = 0.1;
my $time_ini    = time;
my $nread       = 0;
my $index_dir   = '/proj/hoodlab/share/programs/blat-indexes'; # Blat indexes dir
my $fasta       = 'read.fa';
my $psl         = 'read.psl';
my $gfclient    = '/proj/hoodlab/share/programs/blat/gfClient';
my $host        = 'localhost';
my @indexes     = ( 'solexa_primers.2bit',
                    'human.GRCh37.61.rRNA-MT.2bit',
                    'human_RepBase15.10.2bit',
                    'ERCC_reference_081215.2bit',
                    'hs37.61.2bit'                     );
my %indexes     = ();
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
	's|minscore:i'     => \$minscore,
	'p|percent:i'      => \$minident,
	'n|nreads:i'       => \$block,
	'm|maxintron:i'    => \$maxintron,
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
		$nread++;
        my ($id, $seq, $sep, $qual) = split (/\n/, $_);
        $id =~ s/\@//;
        #writeFa($id, $seq);
        my $hit = searchHit($id, $seq);
        print "$hit\n";
		if ($nread % $block == 0) {
			my $time = time - $time_ini;
			warn "query $nread reads in $time seconds\n" if (defined $verbose);
		}
    }
    #unlink "$fasta";
    #unlink "$psl";
}
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
		$nread++;
        my ($id, $seq) = split (/\n/, $_);
        $id =~ s/>//;
        #writeFa($id, $seq);
        my $hit = searchHit($id, $seq);
        print "$hit\n"; 
		if ($nread % $block == 0) {
			my $time = time - $time_ini;
			warn "processed $nread reads in $time seconds\n" if (defined $verbose);
		}
    }
}
else {
    die "Format $format is not FASTQ/FASTA\n";
}

my $time_end = time;
my $time_dif = $time_end - $time_ini;
warn "Query $nread reads take $time_dif seconds\n" if (defined $verbose);

# SUBROUTINES

sub searchHit {
    my ($id, $seq) = @_;
    my $res = "$id\t$seq\t-\t0\tNo_hit_found";
    foreach my $target (@indexes) {
        my $name = $indexes{$target}{'name'};
        my $hit = runBlat($target, ">$id\n$seq\n");
        unless ($hit =~ m/No_hit_found/) {
            $res = "$id\t$seq\t$name\t$hit";
            last;
        }
    }
    return $res;
}

sub runBlat {
    my ($target, $fa) = @_;
    my $port     = $indexes{$target}{'port'};
    my $res      = `echo \"$fa\" | $gfclient $host $port / -nohead -minScore=$minscore -minIdentity=$minident -maxIntron=$maxintron stdin stdout`;
    my @hits     = split (/\n/, $res);
    my $best_hit = 'No_hit_found';
	my $best     = -1;
	my $nhit     = 0;
	# Filter hits, keep the best
    foreach my $hit (@hits) {
        next unless ($hit =~ m/\d+/);
	    chomp $hit;
        my @array = split (/\t/, $hit);
        # Score = (2 * Matches) - Mismatches - QueryGaps - TargetGaps
        my $score = (2 * $array[0]) - $array[1] - $array[4] - $array[6];
		if ($score > $best) {
            $best_hit  = join ":", @array;
            $best = $score;
			$nhit = 1;
        }
        elsif ($score == $best) {
            $best_hit .= '|';
            $best_hit .= join ":", @array;
			$nhit++;
        }
		else {
			# do nothing
		}
    }
    return "$nhit\t$best_hit";
}

sub printVersion {
	print "$0 $our_version\n";
	exit 1;
}
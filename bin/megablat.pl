#!/usr/bin/perl

=head1 NAME

megablat.pl

=head1 DESCRIPTION

Align zillions of short reads to a reference using Blat.
This script will process all the sequences (gfClient) in a file post-
processing of each output.

=head1 USAGE

perl megablat.pl [OPTIONS]

OPTIONS
    Parameters         Description                       Value      Default
    -i --input         Input file  [1]                   FILE       STDIN
    -o --output        Output file [2]                   FILE       STDOUT
    -f --format        Input format                      fq/fa      fq
    -n --nreads        Report time every N reads         INT        1000
    -e --execdir       Directory with the exec           PATH       
    -s --score         Minimal Blat score [3]            INT        100
    -p --percent       Minimal identity percent          INT        90
    -m --maxintron     Maximal intron size               INT        750000
    -r --repeat        Don't report hits more than this  INT        5
    -a --allhits       Keep all hits (no filter step)
    -d --random        Report one random hit
    -x --port          Use this port                     INT        1111
    -y --host          Use this host name                NAME       localhost

    -h --help          Print this screen
    -v --verbose       Activate verbose mode
       --version       Print version number
    
[1] Input file can be compressed (.gz|.bz2)
[2] Output is a PSLX format file, if file exists it will search the last read
    and continue the streaming from that point (crash recovery).
[3] Blat score = (2 * NumMatches) - NumMismatches - NumGaps 

=head1 EXAMPLES

perl megablat.pl -i reads.fq.gz -o align.out

perl megablat.pl -i reads.fa -f fa -o align.out -s 50 -p 95 -n 10000

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
my $host      = 'localhost';   # local host name
my $port      = 1111;          # local port number
my $exec      = '/proj/hoodlab/share/programs/blat';
my $allhits   = undef;         # keep all hits flag
my $repeat    = 5;
my $random    = undef;

# Main variables
my $our_version = 0.1;
my $time_ini    = time;
my $nread       = 0;
my $out         = undef;
my $recover     = 0;

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
    'x|port:i'         => \$port,
    'y|host:s'         => \$host,
    'e|execdir:s'      => \$exec,
    'a|allhits'        => \$allhits,
    'r|repeat:i'       => \$repeat,
    'd|random'         => \$random,
    'version'          => \$version
    ) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if(defined $version);

my $gfclient    = "$exec/gfClient";

# Opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
	if (-e $output) { # simple crash recovery mode
		open F, "$output" or die "output file exits, but I cannot check it\n";
		while (<F>) { 
			$recover++; 
		}
		close F;
		open STDOUT, ">>$output" or die "cannot open $output\n";
	} 
	else { 
		open STDOUT, ">$output" or die "cannot open $output\n";
	}
}

# Main program
if ($format eq 'fq') {
    $/ = "\n\@";
    while (<>) {
        $nread++;
		
		if ($recover > 0) {
			next if ($nread <= $recover);
		}
		
        my ($id, $seq, $sep, $qual) = split (/\n/, $_);
        $id   =~ s/\@//;
        my ($nhit, $hit) = runBlat(">$id\n$seq\n");
        $hit  = '-' if ($nhit >= $repeat);
        $out .= "$id\t$seq\t$nhit\t$hit\n";
        if ($nread % $block == 0) {
            print $out;
            $out     = '';
            my $time = time - $time_ini;
            warn "query $nread reads in $time seconds\n" if (defined $verbose);
        }
    }
}
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
        $nread++;
        
		if ($recover > 0) {
			next if ($nread <= $recover);
		}
		
		my ($id, $seq) = split (/\n/, $_);
        $id   =~ s/>//;
        my ($nhit, $hit) = runBlat(">$id\n$seq\n");
        $hit  = '-' if ($nhit >= $repeat);
        $out .= "$id\t$seq\t$nhit\t$hit\n";
        if ($nread % $block == 0) {
            print $out;
            $out     = '';
            my $time = time - $time_ini;
            warn "processed $nread reads in $time seconds\n" if (defined $verbose);
        }
    }
}
else {
    die "Format $format is not FASTQ/FASTA\n";
}

print $out;
my $time_end = time;
my $time_dif = $time_end - $time_ini;
warn "Query $nread reads take $time_dif seconds\n" if (defined $verbose);

# SUBROUTINES

sub runBlat {
    my $fa       = shift @_;
    my $res      = `echo \"$fa\" | $gfclient $host $port / -nohead -out=pslx -minScore=$minscore -minIdentity=$minident -maxIntron=$maxintron stdin stdout`;
    my @hits     = split (/\n/, $res);
    my $best_hit = 'No_hit_found';
    my $best     = -1;
    my $nhit     = 0;
    
    if (defined $allhits) {
        if (defined $hits[0]) { 
            $nhit = $#hits + 1;
	        $best_hit = join "|", @hits;
        }
	    return "$nhit\t$best_hit";
    }
    # Filter hits, keep the best
    foreach my $hit (@hits) {
        next unless ($hit =~ m/\d+/);
        chomp $hit;
        my @array = split (/\t/, $hit);
        # Score = (2 * Matches) - Mismatches - QueryGaps - TargetGaps
        my $score = (2 * $array[0]) - $array[1] - $array[4] - $array[6];
        if ($score > $best) {
            $best_hit  = join ";", @array;
            $best = $score;
            $nhit = 1;
        }
        elsif ($score == $best) {
            $best_hit .= '|';
            $best_hit .= join ";", @array;
            $nhit++;
        }
        else {
            # do nothing
        }
    }
    return ($nhit, $best_hit);
}

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

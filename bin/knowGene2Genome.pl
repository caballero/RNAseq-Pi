#!/usr/bin/perl

=head1 NAME

knowGene2Genome.pl

=head1 DESCRIPTION

Convert the coordinates from transcript-space to genome-space in a SAM output.
If we have more than one match, it tries to collapse the coordinates (several 
positions in different transcripts could be originated from the same genomic 
location). 

Altered fields are: TARGET, POS_TARGET, CIGAR

=head1 USAGE

perl knowGene2Genome.pl -g human.gtf [OPTIONS] < SAM > SAM

OPTIONS:
   Parameter        Description                    Values       Default
   -i  --input      Read sequences from here       File**       STDIN
   -o  --output     Write formated sequences here  File         STDOUT
   -g  --gtf        Read gene information from GTF File         none
   -x  --excluded   Excluded sequences             File         none
   -v  --verbose    Verbose mode                              
   -h  --help       Print this screen

    ** File can be compressed (gzip/bzip2) but requires option -i
    
=head1 EXAMPLES

perl knowGene2Genome.pl -g human.gtf [OPTIONS] < SAM > SAM

perl knowGene2Genome.pl -g human.gtf -i SAM -o SAM

perl knowGene2Genome.pl -g human.gtf -i SAM.gz > SAM

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
my $help     = undef;         # print help
my $verbose  = undef;         # verbose mode
my $input    = undef;         # input SAM
my $output   = undef;         # output SAM
my $gtf_file = undef;         # annotation GTF
my $excluded = undef;         # excluded sequences

# Main variables
%gtf = ();                    # transcript structure information

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output,
    'g|gtf=s'          => \$gtf_file,
    'e|excluded:s'     => \$exlcuded
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $gft_file);

# opening files (if required)
open GTF, "$gtf_file" or die "cannot open $gtf_file\n";

if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}
if (defined $excluded) {
    open BAD, ">$excluded" or die "cannot write file $excluded\n";
}

warn "loading transcript information from $gft_file\n" if (defined $verbose);
while (<GTF>) {
    my @line = split (/\t/, $_);
    next unless ($line[2] eq 'exon');
    next unless (m/transcript_id "(.+?)"/);
    my $tid = $1;
    my $chr = $line[0];
    my $ini = $line[3];
    my $end = $line[4];
    my $dir = $line[6];
    $gtf{$tid}{'chr'}   = $chr;
    $gtf{$tid}{'dir'}   = $dir;
    my @trans = ();
    for (my $i = $ini; $i <= $end; $i++) {
        push @{ $gtf{$tid}{'trs'} }, $p;
    }
}

# ordering bases in transcripts
foreach my $tid (keys %gtf) {
    my @bases = @{ $gtf{$tid}{'trs'} };
    @{ $gtf{$tid}{'trs'} } = sort { $a<=>$b } (@bases);
}

# parsing the SAM file
while (<>) {
    next if (m/^\@/); # skip SAM headers
    my @line = split (/\t/, $_);
    my $read = $line[0];
    my $flag = $line[1];
    my $hit  = $line[2];
    my $pos  = $line[3];
    my $cig  = $line[5];
    my $dir  = '+'; $dir = '-' if ($flag % 16 == 0);
    if (defined $gtf{$hit}{'chr'}) {
        my ($new_hit, $new_pos, $new_cig, $new_dir) = decodeMap($hit, $pos, $cig, $dir);
        next if (defined $redundant{"$read:$new_hit:$new_pos:$new_cig"});
        $redundant{"$read:$new_hit:$new_pos:$new_cig"}++;
        if ($new_dir eq '+') { $line[1] = 0; } else { $line[1] = 16; }
        $line[2] = $new_hit;
        $line[3] = $new_pos;
        $line[5] = $new_cig;
        $_ = join ("\t", @line);
        print $_;
    } else {
        print BAD $_ if (defined $excluded);
    }
}

sub decodeMap {
    my ($hit, $pos, $cig, $dir) = @_;
    my ($nhit, $npos, $ncig, $ndir);
    $nhit = $gtf{$hit}{'chr'};
    
    my $odir = $gtf{$tid}{'dir'};
    if    ($odir eq '+' and $dir eq '+') { $ndir = '+'; }
    elsif ($odir eq '-' and $dir eq '-') { $ndir = '+'; }
    elsif ($odir eq '+' and $dir eq '-') { $ndir = '-'; }
    elsif ($odir eq '-' and $dir eq '+') { $ndir = '-'; }
    else { die "error comparing directions for $hit:$pos:$dir:$cig\n"; }
    
    my @exons = @{ $gtf{$h}{'trs'} };
    $npos = $exons[$pos];
    my @cig   = split (/(\d+\w)/, $cig);
    my $n = $pos;
    foreach my $slice (@cig) {
        next unless ($splice =~ m/\d+\w/);
        if ($slice =~ m/(\d+)M$/) {
                    }
        elsif ($slice =~ m/(\d+)[IDN]$/) {
            $ncig .= $slice;
        }
        else {
            
        }
    }
    return ($nhit, $npos, $ncig, $ndir);
}

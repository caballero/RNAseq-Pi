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
my %gtf       = ();           # transcript structure information
my %redundant = ();           # remove redudant hits

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output,
    'g|gtf=s'          => \$gtf_file,
    'e|excluded:s'     => \$excluded
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $gtf_file);

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

warn "loading transcript information from $gtf_file\n" if (defined $verbose);
while (<GTF>) {
    my @line = split (/\t/, $_);
    next if ($line[1] =~ m/psudogene/i);
    next unless ($line[2] eq 'exon');
    next unless (m/transcript_id "(.+?)"/);
    my $tid = $1;
    my $chr = $line[0];
    next unless ($chr =~ m/^chr\d+$/ or $chr =~ m/^chr[MYX]$/); # skip other haplotypes
    my $ini = $line[3];
    my $end = $line[4];
    my $dir = $line[6];
    $gtf{$tid}{'chr'}   = $chr;
    $gtf{$tid}{'dir'}   = $dir;
    my @pos = ();
    push @pos, $ini .. $end;
    $gtf{$tid}{'trs'}  .= join ":", @pos;
    $gtf{$tid}{'trs'}  .= ":";
}

# ordering bases in transcripts
warn "sorting transcripts\n" if (defined $verbose);
foreach my $tid (keys %gtf) {
    $gtf{$tid}{'trs'} =~ s/:$//;
    my @bases = split (/:/, $gtf{$tid}{'trs'});
    if ($gtf{$tid}{'dir'} eq '+') {
        @bases = sort { $a<=>$b } (@bases);
    } else {
        @bases = sort { $b<=>$a } (@bases);
    }
    $gtf{$tid}{'trs'} = join ':', @bases;
}

# parsing the SAM file
warn "parsing SAM file\n" if (defined $verbose);
while (<>) {
    chomp;
    next if (m/^\@/); # skip SAM headers
    my @line = split (/\t/, $_);
    my $read = $line[0];
    my $flag = $line[1];
    my $hit  = $line[2];
    my $pos  = $line[3];
    my $cig  = $line[5];
    next unless ($cig =~ m/^\d+M$/); # only exact matches for now (no indels/masking)
    unless (defined $gtf{$hit}{'trs'}) {
        #warn "undefined exon information for $read $hit $pos $cig\n" if (defined $verbose);
        next;
    }
    my $dir  = '+'; 
    $dir = '-' if ($flag == 16);
    if (defined $gtf{$hit}{'chr'}) {
        my ($new_hit, $new_pos, $new_cig, $new_dir) = decodeMap($hit, $pos, $cig, $dir);
        #next if (defined $redundant{"$read:$new_hit:$new_pos:$new_cig"});
        #$redundant{"$read:$new_hit:$new_pos:$new_cig"}++;
        if ($new_dir eq '+') { 
            $line[1] = 0; 
        } 
        else { 
            $line[1] = 16; 
        }
        $line[2] = $new_hit;
        $line[3] = $new_pos;
        $line[5] = $new_cig;
        push @line, "YT:Z:$hit";
        $_ = join ("\t", @line);
        print "$_\n";
    } else {
        print BAD $_ if (defined $excluded);
    }
}

sub decodeMap {
    my ($hit, $pos, $cig, $dir) = @_;
    my ($nhit, $npos, $ncig, $ndir);
    $nhit = $gtf{$hit}{'chr'};
    
    my $odir = $gtf{$hit}{'dir'};
    if    ($odir eq '+' and $dir eq '+') { $ndir = '+'; }
    elsif ($odir eq '-' and $dir eq '-') { $ndir = '+'; }
    elsif ($odir eq '+' and $dir eq '-') { $ndir = '-'; }
    elsif ($odir eq '-' and $dir eq '+') { $ndir = '-'; }
    #else { die "error comparing directions for $hit:$pos:$dir:$cig\n"; }
    
    my @exons = split (/:/, $gtf{$hit}{'trs'});
    my $len   = $cig; 
    $len      =~ s/M$//;
    my $ini   = $pos;
    my $end   = $ini + $len;
    my @ex    = ();
    for (my $i = $ini; $i <= $end; $i++) { 
        push @ex, $exons[$i];
    }
    my $arr1 = length @exons;
    my $arr2 = length @ex;
    
    @ex   = reverse (@ex) if ($odir eq '-');
    $ini  = $ex[1];
    $end  = $ex[-1];
    $npos = $ini;

    warn "something wrong with ini=$ini end=$end len=$len $hit $pos $cig $dir\n" unless (defined $ini and defined $end and defined $len); 
    if ($end - $ini == $len) {
        $ncig = $cig;
    } else {
        my $m = 0;
        for (my $i = 0; $i <= $#ex - 1; $i++) {
            my $diff = $ex[$i + 1] - $ex[$i];
            if ($diff == 1) {
                $m++;
            } 
            else {
                $ncig .= $m . 'M' . $diff . 'N';
                $m = 0;
            }
        }
        $m++;
        $ncig .= $m . 'M';
    }
    return ($nhit, $npos, $ncig, $ndir);
}

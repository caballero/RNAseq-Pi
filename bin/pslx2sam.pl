#!/usr/bin/perl

=head1 NAME

pslx2sam.pl

=head1 DESCRIPTION

Convert a Blat output (modified PSLX format) into a SAM compatible output.

=head1 USAGE

pslx2sam.pl [OPTIONS]

OPTIONS
    Parameters         Description              Values      Default
    -i --input         Input file*              FILE        STDIN
    -o --output        Output file              FILE        STDOUT
    -u --unmap         Keep unmaped reads                   No
    -m --max-hits      Report up to Max hits    INT         10
    -p --poly-map      Report polymapped reads              No
    -s --score         Quality score            Phred+64    I 
    -q --mapq          Map quality score        Phred       255
    
    -h --help          Print this screen
    -v --verbose       Verbose mode
    --version          Print version
  
    * File can be compressed (gzip/bzip2)
  
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

# Main variables
my $version = 0.01;
my $input_h;
my ($sid, $seq, $nhit, $hits);
my %count   = ();

# Parameters initialization
my $getversion = undef;
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $unmap      = undef;
my $polymap    = undef;
my $qual       = 'I';
my $mapq       = 255;
my $maxhits    = 10;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'u|unmap'         => \$unmap,
    'm|max-hits:i'    => \$maxhits,
    'p|polymap'       => \$polymap,
    's|score:s'       => \$qual,
    'q|mapq:i'        => \$mapq,
    'version'         => \$getversion
);
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $getversion);

if (defined $input) {
    $input_h = $input;
    $input_h = "gunzip  -c $input | " if ($input =~ m/\.gz$/);
    $input_h = "bunzip2 -c $input | " if ($input =~ m/\.bz2$/);
    open STDIN, "$input_h" or die "cannot open $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "cannot open $output\n";
}

while (<>) {
    chomp;
    ($sid, $seq, $nhit, $hits) = split (/\t/, $_);
    $count{'total'}++;
    
    if ($hits =~ m/_hap/) {
        filterHap(\$nhit, \$hits);
    }
    
    if ($nhit == 0) {
        printUnmap(\$sid, \$seq, \$nhit) if (defined $unmap);
        $count{'unmap'}++;
    }
    elsif ($nhit > $maxhits) {
        printUnmap(\$sid, \$seq, \$nhit) if (defined $polymap);
        $count{'polymap'}++;
    }
    else {
        printSAM(\$sid, \$seq, \$nhit, \$hits);
        $count{'map'}++;
    }
    
}

warn "Total reads:      $count{'total'}\n" if (defined $verbose);
warn "Unmapped reads:   $count{'unmap'}\n" if (defined $verbose);
warn "Polymapped reads: $count{'poly'}\n"  if (defined $verbose);
warn "Mapped reads:     $count{'map'}\n"   if (defined $verbose);

#####################################################################
#
#                  S U B R O U T I N E S
#
#####################################################################

sub printVersion {
    print "pslx2sam.pl version $version\n";
    exit 1;
}

sub printUnmap {
    my ($sid_ref, $seq_ref, $nhit_ref) = @_;
    my $q = $qual x length $$seq_ref;
    print join "\t", $$sid_ref,4,'*',0,0,'*','*',0,0,$$seq_ref,$q,"NH:i:$$nhit_ref"; 
    print "\n";
}

sub printSAM {
    my ($sid_ref, $seq_ref, $nhit_ref, $hits_ref) = @_;
    my ($chr, $dir, $pos, $read, $q, $mis, $cig, $ctag, $hit);

    
    my @hits = split (/\|/, $$hits_ref);
    foreach $hit (@hits) {
        my @arr = split (/:/, $hit);
        if ($arr[8] eq '-') {
            $dir  = 16;
            $ctag = "XS:A:-";
            $read = revcomp($$seq_ref) if ($dir == 16);
        }    
        else {
            $dir  = 0;
            $ctag = "XS:A:+";
            $read = $$seq_ref;
        }
        $mis  = $arr[1];
        $chr  = $arr[13];
        $pos  = $arr[15] + 1; # 1-based coordinates
        $cig  = '';
        
        $arr[18] =~ s/,$//;
        $arr[19] =~ s/,$//;
        $arr[20] =~ s/,$//;
        
        ## BUILDING CIGAR
        my @blks = split (/,/, $arr[18]); # block lenghts
        my @qblk = split (/,/, $arr[19]); # query ini position in blocks
        my @tblk = split (/,/, $arr[20]); # target ini position in blocks
        my @rpos = split (//, $read);     # all bases slots in read 
        my @qpos = ();
        my @tpos = ();
        my $tpos = 0;
        my $qpos = 0;
        my $dif  = 0;
        # build alignment coordinates puzzle
        foreach my $blk (@blks) {
            $qpos = shift @qblk;
            $tpos = shift @tblk;
            for (1..$blk) {
                $qpos[$qpos] = $qpos + 1;
                $tpos[$qpos] = $tpos + 1;
                $qpos++;
                $tpos++;
            }
        }
        # decode the arrays to obtain the CIGAR
        for (my $i = 0; $i <= $#rpos; $i++) {
            if (defined $qpos[$i] and defined $tpos[$i]) {
                $rpos[$i] = 'M';
            }
            else {
                $rpos[$i] = 'I';
            }
        }
        # soft clip the first bases
        for (my $i = 0; $i <= $#rpos; $i++) {
            last unless ($rpos[$i] eq 'I');
            $rpos[$i] = 'S';
        }
        # soft clip the last bases
        for (my $i = $#rpos; $i >= 0; $i--) {
            last unless ($rpos[$i] eq 'I');
            $rpos[$i] = 'S';
        }
        # calling deletions/splices in the target  
        $tpos = undef;
        for (my $i = 0; $i <= $#rpos; $i++) {
            next unless ($rpos[$i] eq 'M');
            if (defined $tpos) {
                $dif = $tpos[$i] - $tpos;
                if ($dif == 1) {
                    $tpos = $tpos[$i];
                    next;
                }
                else {
                    if ($dif < 70) { # short deletion
                        $rpos[$i] = $dif . 'DM';
                    }
                    else {           # splice site
                        $rpos[$i] = $dif . 'NM';
                    }
                    $tpos = $tpos[$i];
                    next;
                }
            }
            else {
                $tpos = $tpos[$i] if (defined $tpos[$i]); 
            }
        }
        # reconstruct the CIGAR
        my $cod = undef;
        my $cnt = 0;
        push @rpos, 'X';
        for (my $i = 0; $i <= $#rpos; $i++) {
            if ($i == 0) {
                $cod = $rpos[0];
                $cnt = 1;
            }
            else {
                if ($cod eq $rpos[$i]) {
                    $cnt++;
                }
                elsif ($rpos[$i] =~ m/\d+[DN]M/) {
                    $rpos[$i] =~ s/M$//;
                    $cig .= $rpos[$i];
                    $cnt  = 1;
                    $cod  = 'M';
                }
                else {
                    $cig .= $cnt . $cod;
                    $cnt  = 1;
                    $cod  = $rpos[$i];
                }
            }
        }
        
        $q = $qual x length $read;
        print join "\t", $$sid_ref,$dir,$chr,$pos,$mapq,$cig,'*',0,0,$read,$q,"NH:i:".$$nhit_ref,"NM:i:$mis",$ctag; 
        print "\n";
    }
}

sub filterHap {
    my ($nhit_ref, $hits_ref) = @_;
    my $new_nhit = 0;
    my $new_hits = '';
    my @hits = split (/\|/, $$hits_ref);
    foreach my $hit (@hits) {
        next if ($hit =~ m/chr.+hap/);
        $new_nhit++;
        $new_hits .= "$hit|";
    }
    $new_hits =~ s/\|$//;
    $$nhit_ref = $new_nhit;
    $$hits_ref = $new_hits;
}

sub revcomp {
    my $s = shift @_;
    my $r = reverse $s;
    $r =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}

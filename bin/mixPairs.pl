#!/usr/bin/perl

=head1 NAME

mixPairs.pl

=head1 DESCRIPTION

Combine 2 paired-end sequences inputs in one file.

=head1 USAGE

perl mixPairs.pl -1 pair-1.fq -2 pair-2.fq > mixed.fa

OPTIONS:
   Parameter        Description               Values          Default
   -i --in-format   Input file format         fq/fa           fq
   -o --out-format  Output file format        fq/fa           fa
   -m --mixed       Output file               FILE            STDOUT
   -1 --pair1       Input pair 1*             FILE           
   -2 --pair2       Input pair 2*             FILE
   -t --type        Pair relation**           fr/rf/ff        fr
   -q --qual        Default quality score                     B
   -b --batch       Process pairs in batches  INT             1000000
   -v --verbose     Verbose mode
   -h --help        Print this screen
   
   *  File can be compressed (gzip/bzip2)
   ** Pairing can be fr: for-rev, rf: rev-for, ff: for-for. Rev sequences
      will be reverse/complementary.

=head1 EXAMPLES

perl mixPairs.pl -1 pair-1.fq -2 pair-2.fq -i fq -o fa -m mixed.fa -t fr

perl mixPairs.pl -1 pair-1.fq.gz -2 pair-2.fq.gz -i fq -o fq -m mixed.fq -t ff

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
my $help       = undef;    # print help
my $verbose    = undef;    # verbose mode
my $pair1      = undef;    # sequences pair 1
my $pair2      = undef;    # sequences pair 2
my $output     = undef;    # output file
my $format_in  =  'fq';    # input format (fq/fa)
my $format_out =  'fa';    # output format (fq/fa)
my $def_qual   =   'B';    # default quality score
my $type       =  'fr';    # type of pairing
my $batch      =   1e6;    # batch size
my $unpair     = undef;    # write inpair sequences here

# Main variables
my %pairs      = ();
my $cnt        =  0;
my ($id1, $id2, $seq1, $seq2, $sep, $qual1, $qual2);

# Calling options
GetOptions(
    'h|help'         => \$help,
    'v|verbose'      => \$verbose,
    'i|in-format:s'  => \$format_in,
    'o|out-format:s' => \$format_out,
    'm|mixed:s'      => \$output,
    '1|pair1:s'      => \$pair1,
    '2|pair2:s'      => \$pair2,
    't|type:s'       => \$type,
    'b|batch:i'      => \$batch,
    'u|inpaired:s'   => \$unpair,
    'q|qual:s'       => \$def_qual
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $pair1 and defined $pair2);

# Opening files (if required)
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}
if (defined $unpair) {
    open UNP, ">$unpair" or die "cannot open $unpair\n";
}

if ($format_in eq 'fq') {
    my $line = 0;
    my $blob = undef;
    open P1, "$pair1" or die "cannot open $pair1\n";
    while (<P1>) {
        $line++;
        $blob .= $_;
        if ($line == 4) {
            $cnt++;
            ($id1, $seq1, $sep, $qual1) = split (/\n/, $blob);
            $id1 =~ s/\/\d$//;
            $id1 =~ s/^\@//;
            $pairs{$id1}{'s1'} = $seq1;
            $pairs{$id1}{'q1'} = $qual1;
            if ($cnt >= $batch) {
                readPair2();
                flushSeqs();
                $cnt = 0;
            }
            $line = 0;
            $blob = undef;
        }
    }
    readPair2();
    flushSeqs();
}
elsif ($format_in eq 'fa') {
    $/ = "\n\@";
    open P1, "$pair1" or die "cannot open $pair1\n";
    while (<P1>) {
        $cnt++;
        s/>//g;
        ($id1, $seq1) = split (/\n/, $_);
        $id1 =~ s/\/\d$//;
        $pairs{$id1}{'s1'} = $seq1;
        if ($cnt >= $batch) {
                readPair2();
                flushSeqs();
                $cnt = 0;
        }
    }
    readPair2();
    flushSeqs();
}
else {
    die "input format isn't recognized: $format_in\n";    
}
close P1;

flushUnpair() if (defined $unpair);

# SUBROUTINES
sub rc {
    my $s = @_;
    my $r = reverse $s;
    $r =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}

sub readPair2 {
    if ($format_in eq 'fq') {
        open P2, "$pair2" or die "cannot open $pair2\n";
        my $l = 0;
        my $b = undef;
        while (<P2>) {
            $l++;
            $b .= $_;
            if ($l == 4) {
              ($id2, $seq2, $sep, $qual2) = split (/\n/, $b);
              $id2 =~ s/\/\d$//;
              $id2 =~ s/^\@//;
              next unless (defined $pairs{$id2}{'s1'});
              $pairs{$id2}{'s2'} = $seq2;
              $pairs{$id2}{'q2'} = $qual2;
            }
        }
        close P2;
    } else {
        local $/ = "\n>";
        open P2, "$pair2" or die "cannot open $pair2\n";
        while (<P2>) {
              s/>//g;
              ($id2, $seq2) = split (/\n/, $_);
              $id2 =~ s/\/\d$//;
              next unless (defined $pairs{$id2}{'s1'});
              $pairs{$id2}{'s2'} = $seq2;
        }
        close P2;
    }   
}

sub flushSeqs {
    foreach my $id (keys %pairs) {
        if (defined $pairs{$id}{'s1'} and defined $pairs{$id}{'s2'}) {
            $seq1 = $pairs{$id}{'s1'};
            $seq2 = $pairs{$id}{'s2'};
            $qual1 = $def_qual x (length $seq1);
            $qual2 = $def_qual x (length $seq2);
            $qual1 = $pairs{$id}{'q1'} if (defined $pairs{$id}{'q1'});
            $qual2 = $pairs{$id}{'q2'} if (defined $pairs{$id}{'q2'});

            if ($type eq 'fr') {
                $seq2  = rc($seq2);
                $qual2 = reverse($qual2);
            }
            elsif ($type eq 'rf') {
                $seq1  = rc($seq1);
                $qual1 = reverse($qual1);
            }
            elsif ($type eq 'ff') {
                #do nothing
            }
            else {
                die "pairing type isn't recognized: $type\n";
            }
            
            if ($format_out eq 'fq') {
                print "\@$id\n$seq1$seq2\n+\n$qual1$qual2\n";
            }
            elsif ($format_out eq 'fa') {
                print ">$id\n$seq1$seq2\n";
            }
            else {
                die "output format isn't recognized: $format_out\n";
            }
            delete $pairs{$id};
        }
    }
}

sub flushUnpair {
    foreach my $id (keys %pairs) {
        if (defined $pairs{$id}{'s1'}) {
            print UNP "$id/1\t$seq1";
            print UNP "\t", $pairs{$id}{'q1'} if (defined $pairs{$id}{'q1'});
            print UNP "\n";
        }
        elsif (defined $pairs{$id}{'s2'}) {
            print UNP "$id/2\t$seq2";
            print UNP "\t", $pairs{$id}{'q2'} if (defined $pairs{$id}{'q2'});
            print UNP "\n";
        }
    }
}

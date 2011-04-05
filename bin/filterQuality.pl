#!/usr/bin/perl

=head1 NAME

filterQuality.pl

=head1 DESCRIPTION

Script to take a file/stream with Fastq/Fasta sequences and filter low quality 
reads, we check for sequences with too many N's and/or low quality scores.

=head1 USAGE

perl filterQuality.pl [OPTIONS] < FASTA/FASTQ > OUTPUT

OPTIONS:
   Parameter           Description                  Values       Default
   -f  --format        Sequence format              fa/fq        fq
   -c  --color         Color space?                 Y/N          N
   -fn --filter_N      Filter by counting N's?      Y/N          Y
   -fq --filter_qual   Filter by average quality?   Y/N          N
   -fa --filter_array  Filter with array scores?    Y/N          N
   -w  --window        Window size                  int          25
   -n  --max_n         Maximal number of N's        int/float*   3
   -q  --min_qual      Minimal average quality      float        4
   -a  --array_qual    Array of qualities           string       BBBBB
   -t  --trim          Remove bad regions in 3'
   -m  --minsize       Minimal size after trimming  int          25
   -i  --input         Read sequences from here     File**       STDIN
   -o  --out           Write good sequence here     File         STDOUT
   -b  --bad           Write bad sequence here      File
   -g  --group         Process sequences in groups  int          100000
   -v  --verbose       Verbose mode                              
   -h  --help          Print this screen
   
   *  If -n|--max_N is > 1, then we use the integer value, if < 1 we consider
      it as a fraction of the sequence, examples:
          -n 10  => means max 10 N's
          -n 0.2 => means max 20% of N's
   ** File can be compressed (gzip/bzip2) but requires option -i

=head1 EXAMPLES

perl filterQuality.pl < FASTQ > GOOD_FASTQ

perl filterQuality.pl -i FASTQ.gz -o GOOD_FASTQ -b BAD_FASTQ

perl filterQuality.pl -fa Y -a DDDDDDDDDDDDDDDDDDDDD -fq N < FASTQ > GOOD_FASTQ

perl filterQuality.pl -v -f fa -n 8 -t -m 20 -w 20 -i FASTA -o GOOD_FASTA

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2010

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
my $format   = 'fq';          # sequence format
my $color    = undef;         # color-space sequences
my $fil_n    = 'Y';           # filter by number of N's
my $fil_q    = 'N';           # filter by minimal average quality score
my $fil_v    = 'N';           # filter by vector of minimal quality scores
my $max_n    = 3;             # Maximal number of N's
my $min_q    = 4;             # Minimal  average quality score
my $wr_bad   = undef;         # Write bad sequence to this file
my $wr_ok    = undef;         # Write good sequence to this file
my $input    = undef;         # Read sequences from here
my $trim     = undef;         # do trimming in 3'
my $minsize  = 25;            # Minimal size after trimming
my $win      = 25;            # Window to use
my $vec_q    = 'B' x 100;     # Vector of minimal quality scores
my $batch    = 100000;        # Sequences by batch
my $version  = undef;         # Version call flag

# Main variables
my $our_version = 0.1;        # Script version number
my $id       = undef;         # Sequence ID
my $seq      = undef;         # Sequence
my $sep      = undef;         # Fastq separator '+'
my $qual     = undef;         # Fastq qualities
my $tseq     = undef;         # Trimmed sequence
my $tqual    = undef;         # Trimmed qualities
my $sel      = undef;         # Flag for bad sequences
my $selN     = undef;         # Flag for bad sequences (number of N's)
my $selQ     = undef;         # Flag for bad sequences (average quality)
my $selV     = undef;         # Flag for bad sequences (vector of qualities)
my @val      = ();            # General score array
my @valN     = ();            # Score array (number of N's)
my @valQ     = ();            # Score array (average quality)
my @valV     = ();            # Score array (vector of qualities)
my %seqs     = ();            # Hash to keep sequences in memory
$seqs{'num'} = 0;             # Sequences counter
my $nbatch   = 0;             # Counter for batches
my $cnt_ok   = 0;             # Counter for good sequences
my $cnt_bad  = 0;             # Counter for bad sequences
my $cnt_trim = 0;             # Counter for trimmed sequences

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'f|format:s'       => \$format,
    'c|color:s'        => \$color,
    'fn|filter_N:s'    => \$fil_n,
    'fq|filter_qual:s' => \$fil_q,
    'fa|filter_aqual:s'=> \$fil_v,
    'n|num_n:f'        => \$max_n,
    'q|min_qual:f'     => \$min_q,
    'a|array_qual:s'   => \$vec_q,
    'o|out:s'          => \$wr_ok,
    'b|bad:s'          => \$wr_bad,
    't|trim'           => \$trim,
    'w|window:i'       => \$win,
    'g|group:i'        => \$batch,
    'm|minsize:i'      => \$minsize,
    'version'          => \$version
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if(defined $version);

die "You didn't select a filter!\n" unless ( $fil_n ne 'Y' or $fil_q ne 'Y' or $fil_v ne 'Y');

$max_n = int($max_n) if ($max_n >= 1); # remove float

# opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $wr_ok) {
    open STDOUT, ">$wr_ok" or die "cannot write file $wr_ok\n";
}
if (defined $wr_bad) {
    open BAD, ">$wr_bad" or die "cannot write file $wr_bad\n";
}

# Create array of minimal scores (if required)
my @vec_q = qual2num($vec_q) if($fil_v eq 'Y');

# Choosing FASTQ as format
if ($format eq 'fq') {
    my $line = 0;
    my $blob = undef;
    while (<>) {
        $line++;
        $blob .= $_;
        if ($line == 4) {
            $blob =~ s/^\@//;
            $seqs{'num'}++;
            ($id, $seq, $sep, $qual) = split (/\n/, $blob);
            
            # Check for valid sequences
            #my $valS = validateSeq($seq);
            #die "Sequence is not valid: $id + $seq\n" if ($valS == 0);
            #my $valQ = validateQual($qual);
            #die "Quality is not valid: $id + $qual\n" if ($valQ == 0);
        
            # Check the sequence        
            $sel  = 1;
            $selN = 1;
            $selQ = 1;
            $selV = 1;
            if (defined $trim) {
                @valN = checkN($seq,  $win) if ($fil_n eq 'Y');
                @valQ = checkQ($qual, $win) if ($fil_q eq 'Y');
                @valV = checkV($qual, $win) if ($fil_v eq 'Y');
                @val  = mixValues(\@valN, \@valQ, \@valV);
                $tseq = trimSeq($seq, $win, @val);
                if (defined $tseq) {
                    if (length $tseq < $minsize) { 
                        $sel = 0;
                    } 
                    else {
                        if (length $seq != length $tseq) {
                            $cnt_trim++;
                            $tqual = substr($qual, 0, length $tseq);
                            #warn "$id:$seq:$tseq\n" if (defined $verbose);
                            #warn "$id:$qual:$tqual\n" if (defined $verbose);
                        }
                        else {
                            $tqual = $qual;
                        }
                    }
                }
                else { 
                    $sel = 0; 
                }
            }
            else {
                if ($fil_n eq 'Y') {
                    @valN  = checkN($seq,  length $seq);
                    $selN  = $valN[0];
                }
                if ($fil_q eq 'Y') {
                    @valQ  = checkQ($qual, length $seq);
                    $selQ  = $valQ[0];
                }
                if ($fil_v eq 'Y') {
                    @valV  = checkV($qual, length $seq);
                    $selV  = $valV[0];
                }
                if ($selN == 1 and $selQ == 1 and $selV == 1) { $sel = 1; } else { $sel = 0; }
            
                $tseq  = $seq;
                $tqual = $qual;
            }
        
            # Sequences is OK
            if ($sel == 1) { 
                $cnt_ok++;
                $seqs{'good'} .= "\@$id\n$tseq\n$sep\n$tqual\n";
            }
            # Sequence is bad
            else {
                $cnt_bad++;
                $seqs{'bad'} .= "\@$id\n$seq\n$sep\n$qual\n";
            }
    
            flushSeqs() if ($seqs{'num'} >= $batch);
            $line = 0;
            $blob = undef;
        }
    }
    flushSeqs();
}
# Choosing FASTA as format, filter by number of N's only
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
        s/^>//g;
        ($id, $seq) = split (/\n/, $_);
        $seqs{'num'}++;
        $sel = 1;
        
        # Check for valid sequences
        #my $val = validateSeq($seq);
        #die "Sequence is not valid: $id + $seq\n" if ($val == 0);
        
        # Check the sequence
        if (defined $trim) {
            @val = checkN($seq, length $seq);
            $tseq = trimSeq($seq, $win, @val);
            #warn "$seq $tseq @val\n" if (defined $verbose);
            
            if (defined $tseq) {
                if (length $tseq < $minsize) {
                    $sel = 0;
                } 
                else {
                    $cnt_trim++;
                }
            }
            else { 
                $sel = 0; 
            }
        } 
        else {
            @val  = checkN($seq, length $seq);
            $sel  = $val[0];
            $tseq = $seq;
        }
        
        # Sequence is OK
        if ($sel == 1) {
            $cnt_ok++;
            $seqs{'good'} .= ">$id\n$tseq\n";
        }
        # Sequence is bad
        else {
            $cnt_bad++;
            $seqs{'bad'} .= ">$id\n$tseq\n";
        }
        
        flushSeqs() if ($seqs{'num'} >= $batch);
    }
    flushSeqs();
}
# Bad format
else { die "Sorry, I don't recognize this format: $format\n"; } 

# Print report in verbose mode 
if (defined $verbose) {
    my $total    = $cnt_ok + $cnt_bad;
    my $per_ok   = sprintf( "%.2f", 100 * $cnt_ok   / $total);
    my $per_trim = sprintf( "%.2f", 100 * $cnt_trim / $total);
    my $per_bad  = sprintf( "%.2f", 100 * $cnt_bad  / $total);
    warn <<__INFO__    
# Total sequences   = $total
# Good sequences    = $cnt_ok ($per_ok\%)
# Trimmed sequences = $cnt_trim ($per_trim\%)
# Bad sequences     = $cnt_bad ($per_bad\%)
__INFO__
;
}

###################################
####   S U B R O U T I N E S   ####
###################################

# flushSeqs -> print sequences
# Call: flushSeqs()
# Return: nothing
sub flushSeqs  {
    $nbatch += $seqs{'num'};
    #warn "processed $nbatch seqs\n" if (defined $verbose);
    print $seqs{'good'} if (defined $seqs{'good'});
    print BAD $seqs{'bad'} if (defined $wr_bad and defined $seqs{'bad'});
    $seqs{'num'}  = 0;
    $seqs{'good'} = '';
    $seqs{'bad'}  = '';
}

# mixValues -> combine boolean arrays
# Call: mixValues( \@array1, \@array2, ...)
# Return: @values
sub mixValues {
    my @v = ();
    foreach my $ref_arr (@_) {
        my @a = @$ref_arr;
        next unless (defined $a[0]);
        if (defined $v[0]) {
            for (my $i = 0; $i <= $#v; $i++) {
                $v[$i] = 0 if ($v[$i] == 1 and $a[$i] == 0);
            }
        }
        else {
            @v = @a;
        }
    }
    return @v;
}
#  trimSeq -> Remove bad regions
#  Call: trimSeq($seq, @val) STRING, ARRAY
#  Return: $new_seq STRING
sub trimSeq {
    my $seq = shift @_;
    my $win = shift @_;
    my @val = @_;
    my $new = undef;
    my $len = length $seq;
    my $pos = 0;
    foreach my $val (@val) {
        last if ($val == 0);
        $pos += $win;
    }
    
    $pos = $len if ($pos > $len);
    
    if ($pos == $len) { 
        return $seq; 
    }
    elsif ($pos == 0) {
        return undef;
    }
    else {
        return substr($seq, 0, $pos);
    }
}

# checkN => Check for number of N's
# call: checkN($seq, $win)
# return: @values
sub checkN {
    my $s    = shift @_;
    my $w    = shift @_;
    my $n    = 0;
    my @v    = ();
    my $max  = $max_n;
    $max     = int($max_n * length $s) if ($max_n < 1); # max_n is a fraction
    for (my $i = 0; $i < length $s; $i += $w) {
        my $f = substr($s, $i, $w);
        if (defined $color) { $n = $f =~ tr/\./N/; } else { $n = $f =~ tr/N/N/; }
        if ($n <= $max) { push @v, 1; } else { push @v, 0; }
    }
    return @v;
}

# checkQ => Check for average quality score
# call: checkQ($qual, $win)
# return: @values
sub checkQ {
    my $q = shift @_;
    my $w = shift @_;
    my @v = ();
    for (my $i = 0; $i < length $q; $i += $w) {
        my $f = substr($q, $i, $w);
        my @q = qual2num($f);
        my $m = calcMean(@q);
        if ($m >= $min_q) { push @v, 1; } else { push @v, 0; }
    }
    return @v;
}

# checkV => Check for quality scores
# call: checkV($qual, $win)
# return: @values
sub checkV {
    my $q    = shift @_;
    my $w    = shift @_;
    my @v    = ();
    my $max  = $max_n;
    $max     = int($max_n * length $q) if ($max_n < 1); # max_n is a fraction
    for (my $i = 0; $i < length $q; $i += $w) {
        my $n = 0;
        my @q = qual2num( substr($q, $i, $w) );
        for (my $j = 0; $j <= length @q; $j++) {
            $n++ if ($q[$j] <= $vec_q[$j]);
        }
        if ($n >= $max) { push @v, 0; } else { push @v, 1; }
    }
    return @v;
}

# qual2num => Convert ASCII quality scores in decimal values (Phred+64 encoding)
# call: qual2num($qual)
# return: @scores
sub qual2num {
    my $q = shift @_;
    my @q = split (//, $q);
    my @s = ();
    foreach my $x (@q) { push @s, ord($x) - 64; }
    return @s;
}

# calcMean -> calculate the arithmetic mean
# call: calcMean(@values)
# return: $mean
sub calcMean {
    my $sum = 0;
    my $num = 0;
    foreach my $x (@_) { 
        $sum += $x;
        $num++;
    }
    return $sum / $num;
}

# validateSeq => check if the sequence is valid DNA [ACGTN or 0123.]
# call: validateSeq($seq)
# return: 1 => ok, 0 => not ok
sub validateSeq {
    my $s  = shift @_;
    my $ok = 1;
    if (defined $color) { $ok = 0 if ($s =~ m/[^0123\.]/); }
    else                { $ok = 0 if ($s =~ m/[^ACGTN]/i); }
    return $ok;
}

# validateQual => check if the scores are valid (Phred+64 encoding)
# call: validateQual($qual)
# return: 1 => ok, 0 => not ok
sub validateQual {
    my $q  = shift @_;
    my $ok = 1;
    $ok = 0 if ($q =~ m/[^\@ABCDEFGHIJKLMNOPQRSTUVWXYZ\[\\\]\^\_\`abcdefgh]/);
    return $ok;
}

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}


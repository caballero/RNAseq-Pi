#!/usr/bin/perl -w

=head1 NAME

filterLowComplexity.pl

=head1 DESCRIPTION

Script to take a file or stream with FastQ/Fasta sequences and filter low
complexity sequences.

=head1 USAGE

perl filterQuality.pl [OPTIONS] < FASTA/FASTQ > OUTPUT

OPTIONS:
   Parameter         Description                    Values       Default
   -f  --format      Sequence format                fa/fq        fq
   -m  --method      Method to compute complexity   *see below   ct
   -l  --limit       Minimal value for complexity   float        0.02
   -i  --input       Input file                     file**       STDIN
   -o  --output      Output file                    file         STDOUT
   -b  --bad         Bad sequences file             file
   -k  --kmer        K-mer size (if required)       int          3
   -w  --window      Window size                    int          25
   -t  --trim        Trim sequence in 3'
   -n  --minsize     Minimal size after trim        int          25
   -g  --group       Process sequences in groups    int          100000
   -v  --verbose     Verbose mode                                N
   -h  --help        Print this screen
      
   *Methods included:
   Symbol  Description                                    K-mer 
   cwf     Wootton & Federhen complexity                  no
   ce      Entropy of the symbols                         no
   ct      Trifonov's complexity                          1-6
   cl      Linguistic complexity                          1-6
   cm      Markov model of kth-order complexity           1-6
   cz      String compression (Zlib) as complexity        no
   
   ** File can be compressed (gzip/bzip2) but requires option -i

=head1 EXAMPLES

  perl filterLowComplexity.pl < FASTQ > FILTERED.FASTQ
  perl filterLowComplexity.pl -f fa FASTA -o FILTERED.FASTQ -b BAD.FASTQ
  perl filterLowComplexity.pl -f fa FASTA -o FILTERED.FASTQ -m cl -k 2 -t

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
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help     = undef;         # Print help
my $verbose  = undef;         # Verbose mode
my $format   = 'fq';          # Sequence format
my $wr_bad   = undef;         # Write bad sequence to this file
my $wr_ok    = undef;         # Write good sequence to this file
my $input    = undef;         # Read sequences from here
my $method   = 'ct';          # Method to compute complexity
my $kmer     = 3;             # Kmer for method
my $win      = 25;            # Window size
my $minsize  = 25;            # Minimal sequence size after trim
my $trim     = undef;         # Trim the sequences in 3' region
my $lim      = 0.02;          # Complexity cut-off
my $batch    = 100000;        # Write to disk in batches
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
my @val      = ();            # Score array
my $k        = 4;             # Alphabet size [ACGT]
my @alphabet = qw/A C G T/;   # Alphabet
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
    'o|out:s'          => \$wr_ok,
    'b|bad:s'          => \$wr_bad,
    'k|kmer:i'         => \$kmer,
    'l|limit:f'        => \$lim,
    'm|method:s'       => \$method,
    'w|window:i'       => \$win,
    't|trim'           => \$trim,
    'n|minsize:i'      => \$minsize,
    'g|group:i'        => \$batch,
    'version'          => \$version
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $version);
pod2usage(-verbose => 2) unless ($method =~ m/ce|cwf|cl|cz|cm|ct/);

# Aditional module for cz
if ($method eq 'cz') { 
#    use Compress::Zlib;
    die "sorry, Zlib is deactivated for now, plase select another method\n";
}

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

# FastQ format
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
            $sel   = 1;
            @val   = ();
            $tseq  = undef;
            $tqual = undef;
        
            # Trim action
            if (defined $trim) {
                @val   = checkComplex($method, $seq, $kmer, $win);
                $tseq  = trimSeq($seq, $win, @val);
                #warn "$seq $tseq @val\n" if (defined $verbose);
            
                if (defined $tseq) {
                    if (length $tseq < $minsize) { 
                        $sel = 0;
                    } 
                    else {
                        $cnt_trim++;
                        $tqual = substr($qual, 0, length $tseq);
                    }
                }
                else { 
                    $sel = 0; 
                }
            }
            else {
                @val   = checkComplex($method, $seq, $kmer, length $seq);
                $sel   = 0 if ($val[0] < $lim);
                $tseq  = $seq;
                $tqual = $qual
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
            $blob = undef;
            $line = 0;
        }
    }
    flushSeqs();
}
# Fasta format
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
        s/^>//;
        $seqs{'num'}++;
        ($id, $seq) = split (/\n/, $_);
        $sel  = 1;
        @val  = ();
        $tseq = undef;
        
        # Trim action
        if (defined $trim) {
            @val = checkComplex($method, $seq, $kmer, $win);
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
            @val  = checkComplex($method, $seq, $kmer, length $seq);
            $sel  = 0 if ($val[0] < $lim);
            $tseq = $seq;
        }
        
        # Sequences is OK
        if ($sel == 1) { 
            $cnt_ok++;
            $seqs{'good'} .= ">$id\n$tseq\n";
        }
        # Sequence is bad
        else {
            $cnt_bad++;
            $seqs{'bad'} .= ">$id\n$seq\n";
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

#  checkComplex -> Function to apply the correct method
#  Call: checkComplex($method, $seq, $kmer, $win)
#  Return: @res ARRAY
sub checkComplex {
    my ($met, $seq, $kmer, $win) = @_;
    my @res = ();
    if    ($met eq  'ce') { @res =  ce(\$seq, \$win);          }
    elsif ($met eq 'cwf') { @res = cwf(\$seq, \$win);          }
    elsif ($met eq  'cz') { @res =  cz(\$seq, \$win);          }
    elsif ($met eq  'cl') { @res =  cl(\$seq, \$win, \$kmer);  }
    elsif ($met eq  'ct') { @res =  ct(\$seq, \$win, \$kmer);  }
    elsif ($met eq  'cm') { @res =  cm(\$seq, \$win, \$kmer);  }
    else { die "Sorry, I don't recognize this method: $met\n"; }
    return @res;
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
        last if ($val < $lim);
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

#  cwf -> Function to calculate the Complexity by Wootton & Federhen values.
#  Call: cwf( \$seq, \$win, ) STRING, NUMBER
#  Return: @values ARRAY
sub cwf {
	my $seq    = shift;
	my $win    = shift;
	my $len    = length $$seq;
	my @values = ();
	my $up = 0;
	for (my $i = 1; $i <= $$win; $i++) { $up += log_k($k, $$win); }
	for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
		my $str = substr ($$seq, $p, $$win);
		my %elm = countWords($str, 1);
		my $r   = 0;
		my $dw  = 0;
		my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

		foreach my $b (keys %elm) {
			next unless ($elm{$b} > 0);
			$dw += log_k($k, $elm{$b}); 
		}
		$r = ($up - $dw) / $tot if($tot > 1);;
		push @values, $r; 
	}
	return @values;
}

#  ce -> Function to calculate the Complexity Entropy values.
#  Call: ce( \$seq, \$win ) STRING, NUMBER
#  Return: @values ARRAY
sub ce {
	my $seq    = shift;
	my $win    = shift;
	my $len    = length $$seq;
	my @values = ();
	for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
		my $str = substr ($$seq, $p, $$win);
		my %elm = countWords($str, 1);
		my $ce  = 0;		
		my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

		foreach my $b (keys %elm) {
			next unless ($elm{$b} > 0);
			my $r = 0; 
			$r    = $elm{$b} / $tot if($tot > 1); 
			$ce  -= $r * log_k($k, $r); 
		}
		push @values, $ce;
	}
	return @values
}

#  cm -> Function to calculate the Complexity in Markov model values.
#  Call: cm( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER
#  Return: @values ARRAY
sub cm {
	my $seq    = shift;
	my $win    = shift;
	my $word   = shift;
	my $len    = length $$seq;
	my @values = ();
	my $m = pot($k, $$word);
	for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
		my $str = substr ($$seq, $p, $$win);
		my %elm = countWords($str, $$word);
		my $cm  = 0;
		my $dw  = $$win - $$word - 1;
		foreach my $b (keys %elm) {
			next unless ($elm{$b} > 0);
			my $r = $elm{$b} / $dw; 
			$cm  -= $r * log_k($m, $r); 
		}
		push @values, $cm;
	}
	return @values;
}

#  cl -> Function to calculate the Complexity Linguistic values.
#  Call: cl( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER
#  Return: @values ARRAY
sub cl {
	my $seq    = shift;
	my $win    = shift;
	my $word   = shift;
	my $len    = length $$seq;
	my @values = ();
	for (my $p = 0; $p <= ($len - $$win); $p += $$win) {
		my $str     = substr ($$seq, $p, $$win);
		my $sum_vl  = 0;
		my $sum_vm  = 0;
		for (my $l = 1; $l <= $$word; $l++) {
			my $vl  = 0;
			my $vm  = 0;
			my $pot = pot($k, $l); 
			if ($pot < ($$win - $l + 1)) { $vm = $pot;           }
			else                         { $vm = $$win - $l + 1; }
			$sum_vm += $vm;
			my %elm = countWords($str, $l);
			foreach my $b (keys %elm) {
				next unless ($elm{$b} > 0);
				$vl++; 
			}
			$sum_vl += $vl;
		}
		my $r = 0;
		$r    = $sum_vl / $sum_vm if ($sum_vm > 0);
		push @values, $r;
	}
	return @values;
}

#  ct -> Function to calculate the Complexity by Trifonov values.
#  Call: ct( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER
#  Return: @values ARRAY
sub ct {
	my $seq    = shift;
	my $win    = shift;
	my $word   = shift;
	my $len    = length $$seq;
	my @values = ();
	for (my $p = 0; $p <= ($len - $$win); $p  += $$win) {
		my $ct  = 1;
		my $str = substr ($$seq, $p, $$win);
		for (my $l = 1; $l <= $$word; $l++) {
			my $vl  = 0;
			my $vm  = 0;
			my $pot = pot($k, $l);
			if ($pot < ( $$win - $l + 1 )) { $vm = $pot;           }
			else                           { $vm = $$win - $l + 1; }
			my %elm = countWords($str, $l);
			foreach my $b (keys %elm) {
				next unless ($elm{$b} > 0);
				$vl++; 
			}
			$ct *= $vl / $vm;
		}
		push @values, $ct;
	}
	return @values;
}

#  cz -> Function to calculate the Compression factor.
#  Call: cz( \$seq, \$win ) STRING, NUMBER
#  Return: @values ARRAY
#sub cz {
#	my $seq    = shift;
#	my $win    = shift;
#	my $len    = length $$seq;
#	my @values = ();
#	for (my $p = 0; $p <= ($len - $$win); $p++) {
#		my $str = substr ($$seq, $p, $$win);
#		my $r   = 0;
#		my $z   = Compress::Zlib::memGzip($str);
#		my $bz  = length $z;
#		$r      = $bz / $$win; 
#		push @values, $r;
#	}
#	return @values;
#}

#  pot -> Function for calculate the exponential of a number.
#  Call: pot( $num, $exp ) NUMBER, NUMBER
#  Return: $res NUMBER
sub pot {
	my $num = shift @_;
	my $exp = shift @_;
	if ($num == 4) {
		if    ($exp == 1) { return       4; }
		elsif ($exp == 2) { return      16; }
		elsif ($exp == 3) { return      64; }
		elsif ($exp == 4) { return     256; }
		elsif ($exp == 5) { return    1024; }
		elsif ($exp == 6) { return    4096; }
		elsif ($exp == 7) { return   16384; }
		elsif ($exp == 8) { return   65536; }
		elsif ($exp == 9) { return  262144; }
		elsif ($exp ==10) { return 1048576; }
	}
	my $res = $num;
	for (my $i = 2; $i <= $exp; $i++) { $res *= $num; }
	return $res;
}

#  log_k -> Function for calculate the logarithm of any base.
#  Call: log_k( $base, $num ) NUMBER, NUMBER
#  Return: $res NUMBER
sub log_k {
	my $base = shift @_;
	my $num  = shift @_;
	my $res  = 0;
	if ($num > 0) { 
		$res = log $num / log $base; 
	}
	else {
		die "Cannot calculate log_k($base, $num)\n";
	}
	return $res;
}

#  countWords -> Function for count words in a sequence.
#  Call: countWords( $seq, $word ) STRING, NUMBER
#  Return: %results HASH (KEYS are the elements)
sub countWords {
	my $seq   = shift @_;
	my $word  = shift @_;
	my $len   = length $seq; 
	my %count = ();
	
	# Init state when word == 1
	foreach my $b (@alphabet) {
		$count{$b} = 0;
	}
	
	for (my $i = 0; $i <= ($len - $word); $i++) {
		my $elem = substr ($seq, $i, $word);
		next if($elem =~ /[^ACGT]/);
		$count{$elem}++;
	}
	return %count;
}

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}


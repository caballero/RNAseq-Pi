#!/usr/bin/perl

=head1 NAME

knowGene2Genome.pl

=head1 DESCRIPTION

Convert the coordinates from transcript-space to genome-space in a SAM output.
If we have more than one match, it tries to collapse the coordinates (several 
positions in different transcripts could be originated from the same genomic 
location). 

Altered fields are: BITWISE_FLAG, TARGET, POS_TARGET, CIGAR

=head1 USAGE

perl knowGene2Genome.pl -g human.gtf [OPTIONS] < SAM > SAM

OPTIONS:
   Parameter         Description                      Values       Default
   -i  --input       Read sequences from here         File*        STDIN
   -o  --output      Write formated sequences here    File         STDOUT
   -g  --gtf         Read gene information from GTF   File         none
   -e  --excluded    Excluded sequences               File         none
   -u  --uniq        Remove redundant mapped reads               
   -t  --target      Add the target name as "YT:"                  
   -s  --sequences   Read preformated sequences       File         
   -w  --write_seqs  Write preformated sequences      File**
   -n  --no-sam      Don't parse SAM file (for -w)
   -v  --verbose     Verbose mode                              
   -h  --help        Print this screen

   *  File can be compressed (gzip/bzip2) but requires option -i
   ** Options and file names for -g and -s are required
    
=head1 EXAMPLES

perl knowGene2Genome.pl -g human.gtf [OPTIONS] < SAM_cDNA > SAM.chr

perl knowGene2Genome.pl -g human.gtf -i SAM.cDNA -o SAM.chr

perl knowGene2Genome.pl -g human.gtf -i SAM_cDNA.gz -u -t > SAM.chr

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
my $help      = undef;        # print help
my $verbose   = undef;        # verbose mode
my $input     = undef;        # input SAM
my $output    = undef;        # output SAM
my $gtf_file  = undef;        # annotation GTF
my $uniq      = undef;        # print uniquelly mapped reads
my $target    = undef;        # add name of the original target name
my $seq_file  = undef;        # pre-formated sequences file
my $excluded  = undef;        # print excluded sequences here
my $no_sam    = undef;
my $wr_seqs   = undef; 

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
    'e|excluded:s'     => \$excluded,
	'u|uniq'           => \$uniq,
	't|target'         => \$target,
	's|sequences:s'    => \$seq_file,
	'n|no-sam'         => \$no_sam,
	'w|write_seqs'     => \$wr_seqs
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $gtf_file or defined $seq_file);

# opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
open STDOUT, ">$output" or die "cannot write file $output\n" if (defined $output);
open BAD, ">$excluded" or die "cannot write file $excluded\n" if (defined $excluded);

if    (defined $gtf_file) { load_gtf(); }
elsif (defined $seq_file) { load_num_seq(); }
else { die "you didn't load the sequence coordinates!\n"; }

exit 1 if (defined $no_sam);
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
	# filtering of the hits
    unless ($cig =~ m/^\d+M$/) { # only exact matches for now (no indels/masking)
		print BAD $_ if (defined $excluded);
		next;
	}
	if ($flag & 4) { # skip unmapped reads
		print BAD $_ if (defined $excluded);
		next;
	}
	
	if ($flag & 1) {                 # if the read is paired
		unless ($flag & 2) {         # skip unless both pair are properly paired
			print BAD $_ if (defined $excluded);
			next;
		}
	}
	
	unless (defined $gtf{$hit}{'trs'}) {
        #warn "undefined exon information for $read $hit $pos $cig\n" if (defined $verbose);
        print BAD $_ if (defined $excluded);
		next;
    }
    my $dir  = '+'; 
    $dir = '-' if ($flag & 16);
    if (defined $gtf{$hit}{'chr'}) {
        my ($new_hit, $new_pos, $new_cig, $new_dir) = decodeMap($hit, $pos, $cig, $dir);
        if (defined $uniq) {
			next if (defined $redundant{"$read:$new_hit:$new_pos:$new_cig"});
			$redundant{"$read:$new_hit:$new_pos:$new_cig"} = 1;
		}
		
		if ($dir ne $new_dir) { # swap positions in reverse hits
			$line[1]  = flipBitDir($line[1]);
			$line[9]  = rc($line[9]);
			$line[10] = reverse($line[10]);
		}
		
        $line[2] = $new_hit;
        $line[3] = $new_pos;
        $line[5] = $new_cig;
		
		if (defined $target) { # add target name
			push @line, "YT:Z:$hit";
		}
		
		push @line, "XS:A:$new_dir"; # Cufflinks specific field
		
        $_ = join ("\t", @line);
        print "$_\n";
    } 
    else {
        print BAD $_ if (defined $excluded);
    }
}

#-----------------------------------------------#
#           S U B R O U T I N E S               #
#-----------------------------------------------#
sub decodeMap {
    my ($hit, $pos, $cig, $dir) = @_;
    my ($nhit, $npos, $ncig, $ndir);
    $nhit = $gtf{$hit}{'chr'};
    
    my $odir = $gtf{$hit}{'dir'};
    if    ($odir eq '+' and $dir eq '+') { $ndir = '+'; }
    elsif ($odir eq '-' and $dir eq '-') { $ndir = '+'; }
    elsif ($odir eq '+' and $dir eq '-') { $ndir = '-'; }
    elsif ($odir eq '-' and $dir eq '+') { $ndir = '-'; }
    else { die "error comparing directions for $hit:$pos:$dir:$cig\n"; }
    
    my @exons = split (/:/, $gtf{$hit}{'trs'});
    my $len   = $cig; 
    $len      =~ s/M$//;
    my $ini   = $pos - 2;
    my $end   = $ini + $len;
    my @ex    = @exons[$ini .. $end];
    @ex       = reverse (@ex) if ($odir eq '-');
    $ini      = $ex[1];
    $end      = $ex[-1];
    $npos     = $ini;

    if (($end - $ini) == $len) {
        $ncig = $cig;
    } 
    else {
        my $m = 0;
        for (my $i = 0; $i <= $#ex - 1; $i++) {
            my $diff = $ex[$i + 1] - $ex[$i];
            if ($diff == 1) {
                $m++;
            } 
            else {
                $m++;
                $ncig .= $m . 'M' . $diff . 'N';
                $m = 0;
            }
        }
        $ncig .= $m . 'M';
    }
    return ($nhit, $npos, $ncig, $ndir);
}

sub rc {
	my $s = shift @_;
	$s = reverse $s;
	$s =~ tr/acgtACGT/tgcaTGCA/;
	return $s;
}

sub load_gtf {
    my $gtf_file_h = $gtf_file;
    $gtf_file_h = "gunzip  -c $gtf_file | " if ($gtf_file =~ m/\.gz$/);
    $gtf_file_h = "bunzip2 -c $gtf_file | " if ($gtf_file =~ m/\.bz2$/);
    open GTF, "$gtf_file_h" or die "cannot open $gtf_file\n";
    
    if (defined $wr_seqs) {
        if (defined $seq_file) {
            open FT, ">$seq_file" or die "cannot open file $seq_file\n";
        }
        else
        {
            die "please define -s\n";
        }
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
		$gtf{$tid}{'trs'} =~ s/^://;
		my @bases = split (/:/, $gtf{$tid}{'trs'});
		if ($gtf{$tid}{'dir'} eq '+') { @bases = sort { $a<=>$b } (@bases); } 
		else                          { @bases = sort { $b<=>$a } (@bases); }
		$gtf{$tid}{'trs'} = join ':', @bases;
		if (defined $wr_seqs) {
		    print FT join "\t", $tid, $gtf{$tid}{'chr'}, $gtf{$tid}{'dir'}, $gtf{$tid}{'trs'};
		    print FT "\n";
		}
	}
	close FT if (defined $wr_seqs);
	close GTF;
}

sub load_num_seq {
    warn "loading transcript information from $seq_file\n" if (defined $verbose);
    my $seq_file_h = $seq_file;
    $seq_file_h = "gunzip  -c $seq_file | " if ($seq_file =~ m/\.gz$/);
    $seq_file_h = "bunzip2 -c $seq_file | " if ($seq_file =~ m/\.bz2$/);
  	open SEQ, "$seq_file_h" or die "cannot open $seq_file\n";
	while (<SEQ>) {
		chomp;
		my ($tid, $chr, $dir, $pos) = split (/\t/, $_);
		$gtf{$tid}{'chr'} = $chr;
		$gtf{$tid}{'dir'} = $dir;
		$gtf{$tid}{'trs'} = $pos;
	}
	close SEQ;
}

sub flipBitDir {
	my $old = shift @_;
	my $new = $old;
	my @bin = dec2bin($old);
	$bin[4] = swapBit($bin[4]);
	if ($bin[0] == 1) {
		$bin[5] = swapBit($bin[5]);
	}
	$new    = bin2dec(@bin);
	return $new;
}

sub swapBit {
	my $x = shift @_;
	if ($x == 1) {
		return 0;
	}
	else {
		return 1;
	}
}

sub dec2bin {
	my $num = shift @_;
	my @bit = ();
	my @val = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024); #11 bits
	foreach my $val (@val) {
		($num & $val) ? push @bit, 1 : push @bit, 0;
	}
	return @bit;
}
sub bin2dec {
	my $sum = 0;
	my $pot = 0;
	foreach my $bit (@_) {
		$sum += 2**$pot if ($bit == 1);
		$pot++;
	}
	return $sum;
}
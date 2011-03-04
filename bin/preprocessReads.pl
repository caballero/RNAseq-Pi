#!/usr/bin/perl -w

=head1 NAME

preprocessReads.pl 

=head1 DESCRIPTION

Perform basic operations (barcoding, masking, removal) to prepare short reads 
before analysis.

=head1 USAGE

perl preprocessReads.pl [OPTIONS] -i FASTQ/FASTA -o OUTPUT

OPTIONS:
  Parameter           Description                   Values      Default
  -b  --barcode       Separate sequence by barcode  Int    [1]  None
  -m  --mask          Mask this positions (N)       String [2]  None
  -r  --remove        Remove this positions         String [2]  None
  -s  --select        Select this positions         String [2]  None
  -f  --format        Sequence format               fa/fq       fq
  -c  --color         Color-space format                        No
  -i  --input         Read sequences from here      File   [3]  STDIN
  -o  --output        Print good sequence here      File   [4]  STDOUT
  -v  --verbose       Verbose mode                              
  -h  --help          Print this screen  
  
  Notes:
  [1] The value represents the length of the barcode, each barcode produces
      one separate file.
  [2] The value must be 1-base positions to remove separated by ",", ranges 
      can be defined with ":" (i.e. 1:50, 25:, :32).
  [3] Input can be compressed with Gzip (*.gz) or Bzip2 (*.bz2).
  [4] If barcoding, then the name is used as prefix.

=head1 EXAMPLES

   1. Trim sequences to 45 bp
      perl preprocessReads.pl -s 1:45 < FASTQ > TRIMMED
      or
      perl preprocessReads.pl -r 46: < FASTQ > TRIMMED
      
   2. Mask position 24,26 in a Fasta file
      perl preprocessReads.pl -m 24,26 -f fa < FASTA > MASKED
      
   3. Separate reads by barcode (first 6 bases)
      perl preprocessReads.pl -b 6 -i FASTQ -o PREFIX
      
   4. Preprocessing some sequence with barcode (6) and selecting from base 10
      perl preprocessReads.pl -b 6 -s 1:6,10: -i FASTQ -o PREFIX
   
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
my $help     = undef;  # Print help
my $verbose  = undef;  # Verbose mode
my $format   =  'fq';  # Sequence format
my $color    = undef;  # Color-space sequences
my $input    = undef;  # Read sequences from here
my $output   = undef;  # Write ouptut here
my $barcode  = undef;  # Use barcoding
my $select   = undef;  # Select this positions
my $remove   = undef;  # Remove this positions
my $mask     = undef;  # Mask this positions

# Variables
my ($id, $seq, $sep, $qual, $code);
my %openbar = ();      # Hash for barcoding
my $numbar  = 100000;  # Write every $numbar sequences in barcode mode

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output,
    'f|format:s'       => \$format,
    'c|color:s'        => \$color,
    'm|mask:s'         => \$mask,
    'r|remove:s'       => \$remove,
    's|select:s'       => \$select,
    'b|barcode:i'      => \$barcode
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
 
unless(defined $remove or defined $select or defined $barcode or defined $mask) {
    die "You haven't select an operation!\nUse \"perl preprocessReads.pl -h\" for more information\n";
}

if (defined $barcode and !defined $output) {
    if (defined $input) {
        $output = $input;
        $output =~ s/\.(fa|fq|fasta|fastq)//i;
    } else {
        die "You need to assign a prefix for the output with -o\n";
    }
}

# opening files (if required)
if (defined $input) {
    $input = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $input = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
    open STDIN, "$input" or die "cannot read file $input\n";
}
if (defined $output and !defined $barcode) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}

# FASTQ format
if ($format eq 'fq') {
    $/ = "\n\@";
    while (<>) {
        s/^\@//;
        ($id, $seq, $sep, $qual) = split (/\n/, $_);
        
        maskBases(\$mask, \$seq)             if (defined $mask);
        removeBases(\$remove, \$seq, \$qual) if (defined $remove);
        selectBases(\$select, \$seq, \$qual) if (defined $select);
        $seq  =~ s/X//g;
        $qual =~ s/X//g;
        
        if (defined $barcode) {
            $code = substr($seq, 0, $barcode);
            $seq  = substr($seq, $barcode);
            $qual = substr($qual, $barcode);
            $openbar{$code}{'cnt'}++;
            $openbar{$code}{'seq'} .= "\@$id\n$seq\n$sep\n$qual\n";
            
            if ($openbar{$code}{'cnt'} >= $numbar) {
                open  BC, ">>$output-$code.fq" or die "cannot write $output-$code.fq\n";
                print BC  $openbar{$code}{'seq'};
                close BC;
                $openbar{$code}{'seq'} = '';
                $openbar{$code}{'cnt'} = 0;
            }
        }
        else {
            print "\@$id\n$seq\n$sep\n$qual\n"; 
        }

    }
    
    # Flush %openbar
    if (defined $barcode) {
        foreach $code (keys %openbar) {
            open  BC, ">>$output-$code.fq" or die "cannot write $output-$code.fq\n";
            print BC  $openbar{$code}{'seq'};
            close BC;
        }
    }
}
# FASTA format
elsif ($format eq 'fa') {
    $/ = "\n>";
    while (<>) {
        s/^>//;
        ($id, $seq) = split (/\n/, $_);
        
        maskBases(\$mask,     \$seq) if (defined $mask);
        removeBases(\$remove, \$seq) if (defined $remove);
        selectBases(\$select, \$seq) if (defined $select);
        $seq =~ s/X//g;
        
        if (defined $barcode) {
            $code = substr($seq, 0, $barcode);
            $seq  = substr($seq, $barcode);
            $openbar{$code}{'cnt'}++;
            $openbar{$code}{'seq'} .= ">$id\n$seq\n";
            
            if ($openbar{$code}{'cnt'} >= $numbar) {
                open  BC, ">>$output-$code.fa" or die "cannot write $output-$code.fa\n";
                print BC  $openbar{$code}{'seq'};
                close BC;
                $openbar{$code}{'seq'} = '';
                $openbar{$code}{'cnt'} = 0;
            }
        }
        else {
            print ">$id\n$seq\n"; 
        }
    }
    
    # Flush %openbar
    if (defined $barcode) {
        foreach $code (keys %openbar) {
            open  BC, ">>$output-$code.fa" or die "cannot write $output-$code.fa\n";
            print BC  $openbar{$code}{'seq'};
            close BC;
        }
    }
}
# Bad format
else { die "Sorry, I don't recognize this format: $format\n"; } 

###################################
####   S U B R O U T I N E S   ####
###################################

# maskBases -> Change bases to "N" (or "." if color-space)
# Call: maskBases(\$mask, \$seq)
# Return: nothing
sub maskBases {
    my ($ref_range, $ref_seq) = @_;
    my $len = length $$ref_seq;
    my @pos = getPositions($ref_range, $len);
    my $n   = 'N'; $n = '.' if (defined $color);
    foreach my $p (@pos) {
        substr($$ref_seq, $p, 1) = $n;
    }
}

# removeBases -> Delete bases in a sequence array
# Call: maskBases(\$remove, \@seq, \@qual)
# Return: nothing
sub removeBases {
    my ($ref_range, $ref_seq, $ref_qual) = @_;
    my $len = length $$ref_seq;
    my @pos = getPositions($ref_range, $len);
    foreach my $p (@pos) {
        substr($$ref_seq,  $p, 1) = 'X';
        substr($$ref_qual, $p, 1) = 'X' if (defined $$ref_qual);
    }
}

# selectBases -> Delete bases not in a sequence array
# Call: selectBases(\$remove, \@seq, \@qual)
# Return: nothing
sub selectBases {
    my ($ref_range, $ref_seq, $ref_qual) = @_;
    my $len = length $$ref_seq;
    my @pos = getPositions($ref_range, $len);
    my %keep = ();
    foreach my $p (@pos) { 
        $keep{$p} = 1; 
    }
    
    @pos = ();
    for (my $i = 0; $i <= $len; $i++) {
        next if (defined $keep{$i});
        push @pos, $i;
    }
    
    foreach my $p (@pos) {
        substr($$ref_seq,  $p, 1) = 'X';
        substr($$ref_qual, $p, 1) = 'X' if (defined $$ref_qual);
    }
}

# getPositions -> Define the positions wanted and convert into 0-base
# Call: getPositions(\$range, $length)
# Return: @pos
sub getPositions {
    my ($ref_range, $length) = @_;
    my @pos = ();
    my @range = split (/,/, $$ref_range);
    foreach my $r (@range) {
        if ($r =~ m/^\d+$/) {
            $r--;
            push @pos, $r;
        }
        elsif ($r =~ m/^(\d+):(\d+)/) {
            my $ini = $1;
            my $end = $2;
            for (my $i = $ini - 1; $i <= $end - 1; $i++) {
                push @pos, $i;
            }
        }
        elsif ($r =~ m/^(\d+):$/) {
            my $ini = $1;
            my $end = $length;
            for (my $i = $ini - 1; $i <= $end - 1; $i++) {
                push @pos, $i;
            }
        }
        elsif ($r =~ m/^:(\d+)$/) {
            my $ini = 1;
            my $end = $1;
            for (my $i = $ini - 1; $i <= $end - 1; $i++) {
                push @pos, $i;
            }
        }
        else {
            die "Bad range format in: $r\n";
        }
    }
    
    return @pos;
}


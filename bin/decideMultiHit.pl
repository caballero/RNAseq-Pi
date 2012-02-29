#!/usr/bin/perl

=head1 NAME

decideMultiHit.pl

=head1 DESCRIPTION

Read a BAM file and filter reads mapped to unwanted regions (i.e. pseudogenes)
and try to resolve multihit reads prefering gene regions.

=head1 USAGE
  
  perl decideMultiHit.pl -b BAM -g GENES -o OUTPUT

  PARAMETER        DESCRIPTION                VALUE       DEFAULT
  -b --bam         BAM file                   FILE        STDIN
  -g --genes       Gene annotation (BED)      FILE
  -u --unwanted    Bad regions (BED)          FILE
  -o --out         Output file                FILE        STDOUT
  
  -h --help        Print this screen
  -v --verbose     Verbose mode
  

=head1 EXAMPLES

  perl decideMultiHit.pl -b BAM -g GENES -o OUTPUT
  
  perl decideMultiHit.pl -g GENES -o OUTPUT < SAM
  
  samtools view BAM | perl decideMultiHit.pl -g GENES > OUTPUT

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

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
use RSutils;

# Default parameters
my $help       = undef;      # Print help
my $verbose    = undef;      # Verbose mode
my $version    = undef;      # Version call flag
my $bam_file   = undef;      # BAM file
my $genes_file = undef;      # Gene BED file
my $out_file   = undef;      # Output file
my $bad_file   = undef;      # Bad regions BED file
my $bin_size   =   1e5;

# Define where is Samtools, adjust if required
#my $samtools = 'samtools view';   # samtools is in PATH
my $samtools = '/proj/hoodlab/share/programs/samtools/samtools view';

# Main variables
my $our_version = 0.1;     # Script version number
my %genes = ();
my %bads  = ();
my @reads = ();
my $RS    = new RSutils;
my $last  = 'first';

my ($filtered, $chr, $ini, $end, $lab, $bin);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'g|genes:s'        => \$genes_file,
    'b|bam:s'          => \$bam_file,
    'u|unwanted:s'     => \$bad_file,
    'o|out:s'          => \$out
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $genes_file and defined $bad_file);
printVersion() if (defined $version);


warn "loading gene coordinates\n" if (defined $verbose);
open G, "$genes_file" or die "cannot read $genes_file\n";
while (<G>) {
    chomp;
	($chr, $ini, $end, $lab) = split (/\t/, $_);
	$bin = int($ini / $bin_size);
	push @{ $genes{$chr}{$bin} }, "$ini\t$end\t$lab";
}
close G;

warn "loading bad region coordinates\n";
open B, "$bad_file" or die "cannot read $bad_file\n";
while (<B>) {
    chomp;
	($chr, $ini, $end, $lab) = split (/\t/, $_);
	$bin = int($ini / $bin_size);
	push @{ $bads{$chr}{$bin} }, "$ini\t$end\t$lab";
}
close B;

if (defined $bam_file) {
    my $samtools = '/proj/hoodlab/share/programs/samtools/samtools';
    open STDIN, "$samtools view $bam | " or die "cannot open $bam\n";
}

while (<>) {
	chomp;
	my ($read) = split (/t/, $_);
	if ($last eq 'first') {
		$last = $read;
		push @reads, $_;
		next;
	}
	
	if ($last eq $read) {
		push @reads, $_;
		next;
	}
	
	$filtered = filterReads(@reads);
	print "$filtered\n" if (defined $filtered);
	
	$last  = $read;
	@reads = ();
	push @reads, $_;
}

# last read in file
$filtered = filterReads(@reads);
print "$filtered\n" if (defined $filtered);
	
sub filterReads {
	my @filter = ();
	my @genes  = ();
	my ($res, $nh, $cig, $pos, $bad, $gene, $inter, $isBad, $isGene);
	foreach my $read (@_) {
		my @a = split (/\t/, $read);
		$chr  = $a[2];
		$ini  = $a[3];
		$cig  = $a[5];
		$end  = $ini;
		$cig  =~ s/\d+[SID]//g;
		my @cig = split (/[MN]/, $cig);
		foreach my $ext (@cig) { 
			$end += $ext;
		}
		$ini--;
		$pos   = [ "$ini\t$end\tread" ];
		$bin      = int($ini / $bin_size);
		
		# filter bad regions
		if (defined $bads{$chr}{$bin}) {
		    $bad   = \@{ $bads{$chr}{$bin} };
		    $inter = $RS->RSintersection($bad, $pos);
		    $isBad = shift @$inter;
		    next if (defined $isBad);
		}
		
		push @filter, $read;
		
		# filter gene hits
		if (defined $genes{$chr}{$bin}) {
		    $gene   = \@{ $genes{$chr}{$bin} };
		    $inter  = $RS->RSintersection($gene, $pos);
		    $isGene = shift @$inter;
		    push @genes, $r if (defined $isGene);
		}
	}
	if (defined $genes[0]) {
		$res = join "\n", @genes;
		$nh  = $#genes + 1;
		$res =~ s/NH:i:\d+/NH:i:$nh/g;
		return $res;
	} 
	elsif (defined $filter[0]) {
		$res = join "\n", @filter;
		$nh  = $#filter + 1;
		$res =~ s/NH:i:\d+/NH:i:$nh/g;
		return $res;
	}
	else {
		return undef;
	}
}
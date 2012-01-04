#!/usr/bin/perl

use strict;
use warnings;
use RSutils;

my $gene_file = 'genes.block';
my $bad_file  = 'hg19_bad_regions.blocks';

my %genes;
my %bads:
my ($filtered, $chr, $ini, $end, $lab);

warn "loading gene coordinates\n";
open G, "$genes_file" or die "cannot read $genes_file\n";
while (<G>) {
    chomp;
	($chr, $ini, $end, $lab) = split (/\t/, $_);
	push @{ %genes{$chr} }, "$ini\t$end\t$lab";
}
close G;

warn "loading bad region coordinates\n";
open B, "$bad_file" or die "cannot read $bad_file\n";
while (<B>) {
    chomp;
	($chr, $ini, $end, $lab) = split (/\t/, $_);
	push @{ %bads{$chr} }, "$ini\t$end\t$lab";
}
close B;

my $RS = new RSutils;

my $last = 'first';
my @reads;

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
	my @r = @_;
	my @f = ();
	my @g = ();
	my ($res, $nh);
	foreach my $r (@r) {
		my @a = split (/\t/, $r);
		my $chr = $a[2];
		my $ini = $a[3];
		my $cig = $a[5];
		my $end = $ini;
		$cig =~ s/\d+[SID]//g;
		my @cig = split (/[MN]/, $cig);
		foreach my $ext (@cig) { 
			$end += $ext;
		}
		$ini--;
		my $pos   = [ "$ini\t$end\tread" ];
		
		# filter bad regions
		my $bad   = \@{ $bads{$chr} };
		my $inter = $RS->RSintersection($bad, $pos);
		my $isBad = shift @$inter;
		next if (defined $isBad);
		
		push @f, $r;
		
		# filter gene hits
		my $gene   = \@{ $genes{$chr} };
		my $inter  = $RS->RSintersection($gene, $pos);
		my $isGene = shift @$inter;
		push @g, $r if (defined $isGene);
	}
	if (defined $g[0]) {
		$res = join "\n", @g;
		$nh  = $#g + 1;
		$res =~ s/NH:i:\d+/NH:i:$nh/g;
		return $res;
	} 
	elsif (defined $f[0]) {
		$res = join "\n", @f;
		$nh  = $#f + 1;
		$res =~ s/NH:i:\d+/NH:i:$nh/g;
		return $res;
	}
	else {
		return undef;
	}
}
#!/usr/bin/perl

use strict;
use warnings;
use RSutils;

$ARGV[1] or die "use filterBadBlocks.pl BLOCKS CONN\n";

# Infiles
my $block_file = shift @ARGV;
my $conn_file  = shift @ARGV;
my $bad_file   = 'hg19_bad_regions.blocks';

# Filehandlers
my $block_h    = $block_file;
$block_h       = "gunzip  -c $block_file | " if ($block_file =~ m/.gz$/);
$block_h       = "bunzip2 -c $block_file | " if ($block_file =~ m/.bz2$/);
my $conn_h     = $conn_file;
$conn_h        = "gunzip  -c $conn_file | " if ($conn_file =~ m/.gz$/);
$conn_h        = "bunzip2 -c $conn_file | " if ($conn_file =~ m/.bz2$/);

# Outfiles
my $block_out  = $block_file;
$block_out     =~ s/\.gz$//;
$block_out     =~ s/sort/filter/;
my $conn_out   = $conn_file;
$conn_out      =~ s/\.gz$//;
$conn_out      =~ s/sort/filter/;

# Main variables
my $RS         = new RSutils;
my %bad        = ();
my %delBlock   = ();

warn "loading bad region coordinates\n";
open B, "$bad_file" or die "cannot read $bad_file\n";
while (<B>) {
    chomp;
	my ($chr, $ini, $end, $lab) = split (/\t/, $_);
	push @{ $bad{$chr} }, "$ini\t$end\t$lab";
}
close B;

warn "parsing blocks\n";
open BI,    "$block_h" or die "cannot open $block_file\n";
open BO, ">$block_out" or die "cannot open $block_out\n";
while (<BI>) {
	chomp;
	my ($block, $chr, $ini, $end, $ave, $sdt) = split (/\t/, $_);
	if (defined $bad{$chr}[0]) {
		$ini--;
		my $pos   = [ "$ini\t$end\tblock" ];
		my $bad   = \@{ $bad{$chr} };
		my $inter = $RS->RSintersection($bad, $pos);
		my $isBad = shift @$inter;
		if (defined $isBad) {
			$delBlock{$block} = 1;
			next;
		}
	}
	print BO "$_\n";
}
close BI;
close BO;

warn "parsing connections\n";
open CI,    "$conn_h" or die "cannot open $conn_file\n";
open CO, ">$conn_out" or die "cannot open $conn_out\n";
while (<CI>) {
	chomp;
	my ($b1, $b2, $n) = split (/\t/, $_);
	next if ($b1 eq $b2);
	next if (defined $delBlock{$b1});
	next if (defined $delBlock{$b2});
	print CO "$_\n";
}
close CI;
close CO;
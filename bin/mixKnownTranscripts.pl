#!/usr/bin/perl

=head1 NAME

mixKnownTranscripts.pl

=head1 DESCRIPTION

=head1 USAGE


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

my $pair1  = undef;
my $pair2  = undef;
my $fusion = undef;
my $output = undef;
my $unpair = undef;

# Global variables
my %reads  = ();
my ($id, $bit, $hit, $paired);
my ($num_par, $num_orph1, $num_orph2, $num_fus);
my ($tot, $tot1, $tot2);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    '1|pair1:s'        => \$pair1,
    '2|pair2:s'        => \$pair2,
    'f|fusion:s'       => \$fusion,
    'u|unpair:s'       => \$unpair,
    'o|output:s'       => \$output
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if     (defined $help);
pod2usage(-verbose => 2) unless (defined $pair1 and defined $pair2);

# opening files
my $pair1_h = $pair1;
$pair1_h = "gunzip  -c $pair1 | " if ($pair1 =~ m/\.gz$/);
$pair1_h = "bunzip2 -c $pair1 | " if ($pair1 =~ m/\.bz2$/);
open P1, "$pair1_h" or die "cannot open $pair1\n";

my $pair2_h = $pair2;
$pair2_h = "gunzip  -c $pair2 | " if ($pair2 =~ m/\.gz$/);
$pair2_h = "bunzip2 -c $pair2 | " if ($pair2 =~ m/\.bz2$/);
open P2, "$pair2_h" or die "cannot open $pair2\n";

if (defined $fusion) {
	open FS ">$fusion" or die "cannot write $fusion\n";
}

if (defined $unpair) {
	open UP ">$unpair" or die "cannot write $unpair\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}

warn "loading $pair1\n" if (defined $verbose);
while (<P1>) {
	next if (m/^\@/);
	my @line = split (/\t/, $_);
	$id   = $line[0];
	$bit  = $line[1];
	next if ($bit == 4);
	$tot1++ unless (defined $reads{$id}{'1'});
	$id =~ s/#.+$//;
	$hit  = $line[3];
	$reads{$id}{'1'}{$hit} = $bit;
}

warn "loading $pair2\n" if (defined $verbose);
while (<P2>) {
	next if (m/^\@/);
	my @line = split (/\t/, $_);
	$id   = $line[0];
	$bit  = $line[1];
	next if ($bit == 4);
	$tot2++ unless (defined $reads{$id}{'2'});
	$id =~ s/#.+$//;
	$hit  = $line[3];
	$reads{$id}{'2'}{$hit} = $bit;
}

warn "wow, I'm still alive, pairing reads\n" if (defined $verbose);
foreach $id (keys %reads) {
	$tot++;
	if (defined $reads{$id}{'1'}) {
		$paired = 0;
		if (defined $reads{$id}{'2'} {
			foreach $hit (keys %{ $reads{$id}{'1'} }) {
				if (defined $reads{$id}{'2'}{$hit}) {
					$paired = 1;
					$num_par++;
					print "$id\t$hit\n";
				}
			}
			if ($paired == 0) {
				$num_fus++;
				if (defined $fusion) {
					print FS "$id\t"; 
					print FS join (":", keys %{ $reads{$id}{'1'} });
					print FS "\t";
					print FS join (":", keys %{ $reads{$id}{'2'} });
					print FS "\n";
				}
			}
		}
		else {
			$num_orph1++;
			if (defined $unpair) {
				print UP "$id\t", join (":", keys %{ $reads{$id}{'1'} }), "\tNA\n";
			}
		}
	}
	else {
		if (defined $unpair) {
			$num_orph2++;
			print UP "$id\tNA\t", join (":", keys %{ $reads{$id}{'2'} }), "\n";
		}
	}
}

if (defined $verbose) {
	warn "Total reads: $tot\nReads in $pair1: $tot1\nReads in $pair2: $tot2\n";
	warn "Paired reads: $num_par\n";
	warn "Fusion detected: $num_fus\n" if (defined $fusion);
	warn "Orphans in $pair1: $num_orph1\nOrphans in $pair2: $num_orph2\n" if (defined $unpair); 
}
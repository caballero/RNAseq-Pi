#!/usr/bin/perl

=head1 NAME

combineSAM.pl 

=head1 DESCRIPTION

Read a list of SAM files, get the mapped reads and select the best align for each read.

=head1 USAGE

perl combineSAM.pl -i SAM1,SAM2,... -o SAM

OPTIONS:
   Parameter       Description                        Values   Default
   -i --input      SAM files (separated with commas)  Files*      
   -o --output     Write sequences here               File     STDOUT
   -v --verbose    Verbose mode
   -h --help       Print help
   
   *File can be compressed (gzip/bzip2), it requires at least 2 SAM files. 
   
=head1 EXAMPLES

perl combineSAM.pl -i align1.sam,align2.sam,align3.sam -o combined_align.sam

perl combineSAM.pl -v -i align1.sam.gz,align2.sam.gz,align3.sam.gz > combined_align.sam

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
my $sam_files = undef;
my $output    = undef;
my $help      = undef;
my $verbose   = undef;

# Main variables
my @sam_files = ();
my %data      = ();

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'i|input:s'        => \$input,
    'o|output:s'       => \$output
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);

@sam_files = split (/,/, $input);
pod2usage(-verbose => 2) if (length @sam_files < 2);

# opening files (if required)
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write file $output\n";
}

foreach my $sam (@sam_files) {
    my $name = $sam;
	$name =~ s/\.sam.*//;
	$sam = "gunzip  -c $input | " if ($sam =~ m/gz$/);
    $sam = "bunzip2 -c $input | " if ($sam =~ m/bz2$/);
    open SAM, "$sam" or die "cannot read file $sam\n";
	while (<SAM>) {
		next if (m/^\@/);
		my @line = split (/\t/, $_);
		my $rid  = $line[0];
		my $maq  = $line[4];
		m/NM:i:(\d+)/;
		my $mis  = $1;
		if (defined $data{$rid}{'aln'}) {
			my $old_maq = $data{$rid}{'maq'};
			my $old_mis = $data{$rid}{'mis'};
			my $old_aln = $data{$rid}{'aln'};
			
			if ($maq > $old_maq) {
				$data{$rid}{'aln'} = $_;
				$data{$rid}{'maq'} = $maq;
				$data{$rid}{'mis'} = $mis;
			}
			elsif ($mis < $old_mis) {
				$data{$rid}{'aln'} = $_;
				$data{$rid}{'maq'} = $maq;
				$data{$rid}{'mis'} = $mis;
			}
			elsif ($maq == $old_maq and $mis == $old_mis and $_ ne $old_aln) {
				$data{$rid}{'aln'} .= $_;
			}
			else {
				#do nothing, keep last record
			}
		}
		else {
			$data{$rid}{'aln'} = $_;
			$data{$rid}{'maq'} = $maq;
			$data{$rid}{'mis'} = $mis;
		}
	}
	close SAM;
}

foreach my $read (keys %data) {
	print $data{$rid}{'aln'};
}
#!/usr/bin/perl

=head1 NAME

filterTopHatFusionSAM.pl

=head1 DESCRIPTION


=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
    -t --top         Header                     File       sam_header
    -h --help        Print this screen
    -v --verbose     Verbose mode

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

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $header     = 'sam_header.txt';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    's|sam:s'         => \$output,
    't|top:s'         => \$header
    
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    my $input_h = $input;
    $input_h = "gunzip  -c $input | " if ($input =~ m/gz$/i);
    $input_h = "bunzip2 -c $input | " if ($input =~ m/bz2$/i);
    open STDIN, "$input_h" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

open HEADER, ">$header" or die "Cannot write file $header\n";

while (<>) {
    if (m/^\@/) {
        #it's a header line
        print HEADER $_;
    }
    else {
        my @line = split (/\t/, $_);
        next if ($line[1] & 4); # Unmapped flag true
        next if ($line[5] =~ m/F/); # Fussion flag in TopHat-fussion
        
        my $dir = '+';
        $dir = '-' if ($line[1] & 16); # reverse flag true
        
        chomp;
        print "$_\tXS:A:$dir\n";
    }
}

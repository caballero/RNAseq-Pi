#!/usr/bin/perl

=head1 NAME

reportMappedReads.pl

=head1 DESCRIPTION

Check all files in a directory, build a report of sequences in SAM and FASTQ
files.

=head1 USAGE

perl reportMappedReads.pl [OPTIONS]

OPTIONS:
    Parameter       Description                 Value       Default
    -d  --dir       Check files in DIR          PATH        local dir
    -n  --name      Analyze files with NAME     NAME        all files
    -o  --output    Write output here           FILE        STDOUT
    -e  --exclude   Exclude files with NAME     NAME        
    -v  --verbose   Verbose mode
    -h  --help      Print this screen

=head1 EXAMPLES

perl reportMappedReads.pl

perl reportMappedReads.pl -n my_sample -o my_sample.report

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

# Option variables
my $help    = undef;
my $verbose = undef;
my $output  = undef;
my $name    = undef;
my $exclude = undef;
my $dir     = '.';

# Global variables
my ($file, $count);
my $types = 'fq|fastq|fa|fasta|sam';
# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'n|name:s'         => \$name,
    'd|dir:s'          => \$dir,
    'e|exclude:s'      => \$exclude,
    'o|output:s'       => \$output
);

pod2usage(-verbose => 2) if (defined $help);

# opening files
opendir DIR, "$dir" or die "cannot open directory $dir\n";

if (defined $output) {
    open STDOUT, ">$output" or die "cannot open $output\n";
}

# reading files in directory
while ($file = readdir DIR) {
    if (defined $name) {
        next unless ($file =~ m/$name/);
    }
    
    if (defined $exclude) {
        next if ($file =~ m/$exclude/);
    }
    
    next unless ($file =~ m/[$types](.gz|.bz2)*$/i);
    my $type  = defineType($file);
    my $count = countReads($file, $type);
    print "$file\t$count\n"
}


# Subroutines
sub defineType {
    my $f = shift @_;
    my $t = 'NA';
    if    ($f =~ m/[fq|fastq](.gz|.bz2)*$/i) { $t = 'FQ'; }
    elsif ($f =~ m/[fa|fasta](.gz|.bz2)*$/i) { $t = 'FA'; }
    elsif ($f =~ m/sam(.gz|.bz2)*$/i)        { $t = 'SAM'; }
    
    return $t;
}

sub  countReads {
    my ($f, $t) = @_;
    my $c = 0;
    my %s = ();
    
    my $fh = $f;
    $fh = "gunzip -c $f | "  if ($f =~ m/.gz$/);
    $fh = "bunzip2 -c $f | " if ($f =~ m/.bz2$/);
    
    open F, "$fh" or die "cannot open $fh\n";
    while (<F>) {
        if ($t eq 'SAM') {
            next if (m/^\@/);
            my @a = split (/\t/, $_);
            next if ($a[1] == 4);
            next if (defined $s{$a[0]});
            $s{$a[0]} = 1;
            $c++;
        }
        elsif ($t eq 'FA') {
            $c++ if (m/^>/);
        }
        elsif ($t eq 'FQ') {
            $c++ if (m/^\@/);
        }
        else {
            warn "format not recognized for $f => $t\n";
            $c = 'NA';
        }
    }
    close F;

    return $c;
}

#!/usr/bin/perl

=head1 NAME

blatServers.pl

=head1 DESCRIPTION

Launch or stop Blat gfServers. 

=head1 USAGE

perl blatServers.pl [OPTIONS]

OPTIONS
    Parameter          Description   
	-s --start         Start servers 
	-x --stop          Stop servers

	-h --help          Print this screen
	-v --verbose       Activate verbose mode
	--version          Print version number
	
=head1 EXAMPLES

Launch servers:
    perl blatServers.pl -s
	
Stop servers:
    perl blatServers.pl -x

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

# Default parameters
my $help      = undef;         # print help
my $verbose   = undef;         # verbose mode
my $version   = undef;         # version call flag
my $start     = undef;         # start servers flag
my $stop      = undef;         # stop servers flag

# Configuration
my $our_version = 0.1;
my $index_dir = '/proj/hoodlab/share/programs/blat-indexes'; 
my $gfserver  = '/proj/hoodlab/share/programs/blat/gfServer';
my $param     = '-canStop';
my $host      = 'localhost';

my %indexes   = ();
$indexes{'solexa_primers.2bit'}          = 1111;
$indexes{'human.GRCh37.61.rRNA-MT.2bit'} = 1112;
$indexes{'human_RepBase15.10.2bit'}      = 1113;
$indexes{'ERCC_reference_081215.2bit'}   = 1114;
$indexes{'hs37.61.2bit'}                 = 1115;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    's|start'          => \$start,
    'x|stop'           => \$stop,
	'version'          => \$version
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $start or defined $stop);

printVersion() if(defined $version);

if (defined $start) {
    foreach my $idx (keys %indexes) {
        my $port = $indexes{$idx};
        warn "launching server for $idx using host=$host, port=$port\n" if (defined $verbose);
        system ("$gfserver start $host $port $param $index_dir/$idx &");
    }
}
else {
    foreach my $idx (keys %indexes) {
        my $port = $indexes{$idx};
		warn "stoping server for $idx in host=$host, port=$port\n" if (defined $verbose);
        system ("$gfserver stop $host $port");
    }
}

sub printVersion {
	print "$0 $our_version\n";
	exit 1;
}
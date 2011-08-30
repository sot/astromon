#!/usr/bin/env perl

use warnings;
use strict;
use IO::All;
use Ska::Run;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long;

our %opt = (select_name => 'standard_xcorr');
GetOptions(\%opt,
           'select_name=s');

$Ska::Run::LOUD = 1;

our $ASTROMON_DATA  = abs_path(dirname(__FILE__)) . "/data";
our $ASTROMON_WWW   = abs_path(dirname(__FILE__)) . "/data/$opt{select_name}/www" ;

chdir "$ASTROMON_DATA/$opt{select_name}";
io($ASTROMON_WWW)->mkpath;

foreach my $det (qw(ACIS-S ACIS-I HRC-S HRC-I)) {
    run("ps2any -size 475 -rotate 90 " .
	"offsets-${det}.ps $ASTROMON_WWW/offsets-${det}.gif");
    run("ps2any -size 400 -rotate 90 " .
	"offsets-${det}-hist.ps $ASTROMON_WWW/offsets-${det}-hist.gif");
}

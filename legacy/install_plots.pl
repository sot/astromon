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
    run("convert -rotate 90 -quality 95 -density 72x72 -geometry 475 -normalize -sharpen 0 "
        . " -flatten offsets-${det}.ps $ASTROMON_WWW/offsets-${det}.gif");
    run("convert -rotate 90 -quality 95 -density 72x72 -geometry 400 -normalize -sharpen 0 "
        . " -flatten offsets-${det}-hist.ps $ASTROMON_WWW/offsets-${det}-hist.gif");
}

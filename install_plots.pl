#!/usr/bin/env /proj/sot/ska/bin/perl

use warnings;
use strict;
use IO::All;
use Ska::Run;

$Ska::Run::LOUD = 1;

our $ASTROMON_DATA  = "$ENV{SKA}/data/astromon";
our $ASTROMON_WWW   = "$ENV{SKA}/www/ASPECT_PUBLIC/celmon";

chdir $ASTROMON_DATA;
io($ASTROMON_WWW)->mkpath;

foreach my $det (qw(ACIS-S ACIS-I HRC-S HRC-I)) {
    run("ps2any -size 375 -rotate 90 " .
	"offsets-${det}.ps $ASTROMON_WWW/offsets-${det}.gif");
    run("ps2any -size 300 -rotate 90  " .
	"offsets-${det}-hist.ps $ASTROMON_WWW/offsets-${det}-hist.gif");
}

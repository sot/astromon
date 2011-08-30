#!/usr/bin/env perl

use warnings;
use strict;

use IO::All;
use Data::Dumper;

my $table_def < io('astromon_table_defs');
my %table = eval $table_def;

print Dumper \%table;
my %short = ('astromon_obs' => 'o',
	     'astromon_xcorr' => 'xc',
	     'astromon_cat_src' => 'c',
	     'astromon_xray_src' => 'x',
	     );

while (my ($long, $short) = each %short) {
    my @field_specs = @{$table{$long}};
    while (@field_specs) {
	my $field = shift @field_specs;
	next if ref($field);
	my $spec  = shift @field_specs;
	print "   $short.$field as ${short}_${field},\n";
    }
}


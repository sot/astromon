#!/usr/bin/env /proj/sot/ska/bin/perlska

use warnings;
use strict;

use PDL;
use PDL::Graphics::PGPLOT::Window;

our $PLOT_DATA_FILE = "$CELMON_DATA/PLOT_DATA.rdb";

##***************************************************************************
sub make_plots {
##***************************************************************************
    $tmp_data_file = $PLOT_DATA_FILE;
    %det_mapping = ('ACIS-S' => 0,
		    'ACIS-I' => 1,
		    'HRC-S' => 2,
		    'HRC-I' => 3);
    
# Read celestial location database

    print "Making celestial location plots\n" if ($par{loud});
    $rdb = new RDB "${cel_loc_db}" or die "Couldn't open ${cel_loc_db}\n";
    while( $rdb->read( \%data )) {
	unless ($data{status}) { # Any non-blank status indicates a problem with data
	    push @db, {%data};
	    $obsid{$data{obsid}} = 1;
	}
    }
    undef $rdb;

    open DATA, "> $tmp_data_file" or die "Couldn't open $tmp_data_file\n";
    print DATA join("\t", qw(detector tstart dy dz pos_ref version obsid det date_obs)), "\n";
    print DATA join("\t", qw(  N        N    N   N    N      N      N     S       s)), "\n";
    foreach $obsid (keys %obsid) {
	@dbo = grep {$_->{obsid} == $obsid} @db;
	@cat = unique_sources(@dbo);
	foreach $cat (@cat) {
	    next unless ($cat->{field_dr} < $par{max_field_dr});
	    print DATA join("\t", $det_mapping{$cat->{detector}},
				  $cat->{tstart},
				  $cat->{dy},
				  $cat->{dz},
				  catalog_order($cat->{pos_ref}),
				  $cat->{version},
				  $obsid,
				  $cat->{detector},
			          $cat->{date_obs}), "\n";
	}
    }
    close DATA;

    plot_time_history($tmp_data_file);

    system("ps2gif -scale 1.3 $CELMON_DATA/offsets.ps");
    print "Files offsets.ps and offsets.gif created\n" if ($par{loud});
}

##*****************************************************************************
sub plot_time_history {
##*****************************************************************************
    my $file = shift;
    my @dets = qw(ACIS-S ACIS-I HRC-S HRC-I);
    my $d = Ska::HashTable->new($file)->row("pos_ref <= 1");

    my $device = "$CELMON_DATA/offsets.ps/cps";	# color postscript file

    my $win = PDL::Graphics::PGPLOT::Window->new(Device => $device,
					      NXPanel => 1,
					      NYPanel => 4,
					      WindowXSize => 7,
					      WindowYSize => 10.
					     );

    my $t = pdl $d->col("tstart");
    my $yr = $t / 86400. / 365.25 + 1998.0;
    my ($x0, $x1) = minmax($yr);
    $x0 -= 0.1;
    $x1 += 0.1;
    my ($y0, $y1) = (-3, 3);

    for my $det (0 .. 3) {
	$win->env($x0, $x1, $y0, $y1, {Panel => $det+1, charsize=>2.0});
	$win->label_axes("Year", "Offset (arcsec)", "$dets[$det]: DY = red triangle  DZ = blue square",
			 {
			  charsize=>2.0});
	my $d_det = $d->row("detector eq $det");
	my $dy = pdl $d_det->col("dy");
	my $dz = pdl $d_det->col("dz");
	my $x  = pdl($d_det->col("tstart")) / 86400. / 365.25 + 1998.0;
	$win->points($x, $dy, {color => 'red', symbol => 13, charsize=>2.5});
	$win->points($x, $dz, {color => 'blue', symbol => 16, charsize=>2.5});

	map { $win->line([$x0,$x1], [$_, $_], {linestyle => 'dotted'}) } (-2,-1,0,1,2);
    }    
}

##*****************************************************************************
sub unique_sources {
##*****************************************************************************
    my @db = @_;
    my @src = ();
    my $match;

    foreach $db (@db) {
	$match = 0;
	foreach $src (@src) {
	    if (radec_dist($db->{ra}, $db->{dec}, $src->{ra}, $src->{dec})*3600 < $par{det_rad}) {
		$match = 1;
		# Replace the source entry if DB version is "better" than source
		$src = $db if ($db->{version} > $src->{version}
			       || catalog_order($db->{pos_ref}) < catalog_order($src->{pos_ref}));
		last;
	    }
	}
	push @src, $db unless ($match);
    }
    return @src;
}

##*****************************************************************************
sub catalog_order {
##*****************************************************************************
    $_ = shift;
    my $ord = 3;

    $ord = 2 if (/2MASS/ || /USNO/ || /CELMON/ || /SIMBAD_med/);
    $ord = 1 if (/SIMBAD_high/ || /CELMON/);
    $ord = 0 if (/ICRS/ || /Tycho/);

#   $ord = 3 if (/SIMBAD/);

    return $ord;
}


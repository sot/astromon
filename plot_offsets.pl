#!/usr/bin/env perl

use warnings;
use strict;

use DBI;
use IO::All;
use Getopt::Long;
use Ska::Process qw(get_params run_tool get_archive_files);
use Ska::DatabaseUtil qw(:all);
use Ska::RDB qw(write_rdb);
use Ska::HashTable;
use Chandra::Time;
use PDL;
use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::LUT;
use PDL::NiceSlice;
use Data::Dumper;
use List::MoreUtils;
use App::Env qw(CIAO);
use Carp;
use Chandra::Tools::dmcoords;
use File::Basename qw(dirname);
use Cwd qw(abs_path);

our $ASTROMON_DATA  = abs_path(dirname(__FILE__)) . "/data";
our %SIM_z_nom = ('acis-s' => -190.14,
                  'acis-i' => -233.59,
                  'hrc-s'  => 250.47,
                  'hrc-i'  => 126.98);

our %evt2_files; # Event files downloaded to make images.  Ask if want to delete on exit.

our %opt = (select_name => 'standard_xcorr',
            batch => 0,
            pos_ref_lim => 11,  # Include SDSS
            sim_offset => 4.0,  # Max SIM-z offset in mm
            angle_lim   => 240,
            snr_lim     => 4.0,
            all         => 0,
            tstart      => '1998:001:00:00:00',
            tstop       => '2037:001:00:00:00',
           );
our %title = (dy_dz => 'DY = red triangle  DZ = blue square',
              dr    => 'DR = green circle',
             );

GetOptions(\%opt,
           'select_name=s',
           'batch!',
           'sim_offset=f',
           'pos_ref_lim=i',
           'angle_lim=f',
           'snr_lim=f',
           'all!',
           'tstart=s',
           'tstop=s',
          );

# Make directory for plots
our $PLOT_DIR = sprintf("$ASTROMON_DATA/$opt{select_name}");
io($PLOT_DIR)->mkpath unless io($PLOT_DIR)->exists;

# Convert start and stop times to secs
$opt{$_} = Chandra::Time->new($opt{$_})->secs() for (qw(tstart tstop));

my $query = <<END_OF_QUERY ;
SELECT *,
       x.ra AS x_ra, x.dec AS x_dec,
       c.ra AS c_ra, c.dec AS c_dec
  FROM  astromon_xcorr AS xc
    JOIN astromon_cat_src  AS c ON xc.obsid = c.obsid AND xc.c_id = c.id
    JOIN astromon_xray_src AS x ON xc.obsid = x.obsid AND xc.x_id = x.id
    JOIN astromon_obs      AS o ON xc.obsid = o.obsid
    WHERE xc.select_name='$opt{select_name}'
       AND x.r_angle<=$opt{angle_lim}
       AND x.snr>=$opt{snr_lim}
       AND o.tstart >= $opt{tstart}
       AND o.tstart < $opt{tstop}
END_OF_QUERY
    ;

# Get sqlite database handle for astromon database
my $dbh = get_dbh("$ASTROMON_DATA/astromon.db3");

my @query_data = sql_fetchall_array_of_hashref($dbh, $query);

# Take the raw cross-correlation data and massage a bit for plotting
my @plot_data = get_plot_data(@query_data);

# Write to an RDB file for convenience later
my $plot_data_file = "$ASTROMON_DATA/$opt{select_name}/plot.rdb";
write_rdb($plot_data_file,
          \@plot_data,
         );

# Now plot using the RDB file
plot_all_offsets($plot_data_file);

# Close down arc5gl properly
undef $Ska::Process::arc5gl;

# Ask user about deleting left over event files
# First make sure the files are still there:
do { delete $evt2_files{$_} unless -e $evt2_files{$_}} for keys %evt2_files;

if (my @evt2_files = keys %evt2_files) {
    print "Some event files were left over:\n";
    print join("\n", @evt2_files), "\n";
    print "Clean now (y/N)? ";
    chomp (my $a = <STDIN>);
    if ($a eq 'y') {
        print "Deleting files..\n";
        unlink @evt2_files{@evt2_files} ;
    }
}

##****************************************************************************
sub get_dbh {
##****************************************************************************
    ##-- Make database connection
    my $dbfile  = shift;
    my $dbh = DBI->connect("DBI:SQLite:dbname=$dbfile", "", "");
    die "Couldn't connect: $DBI::errstr" unless ($dbh);

    return $dbh;
}

##***************************************************************************
sub get_plot_data {
##***************************************************************************
    my @db = @_;
    my @plot_data;
    my %det_ref = ('acis-s' => 0,
                   'acis-i' => 1,
                   'hrc-s' => 2,
                   'hrc-i' => 3);

    my @obsids = List::MoreUtils::uniq map {$_->{obsid}} @db;
    foreach my $obsid (@obsids) {
        my @dbo = grep {$_->{obsid} == $obsid} @db;
        my @cat = unique_sources(@dbo);
        foreach my $cat (@cat) {
            # Annoying hack because Data::ParseTable fails for fields like BAADE'S WINDOW
            foreach (values %{$cat}) {
                s/['"]//g if defined $_;
            }

            # Ignore sources with a SIM-z offset that is too large (more than 4mm by default)
            next if abs($cat->{sim_z} - $SIM_z_nom{ $cat->{detector} }) > $opt{sim_offset};

            push @plot_data, { pos_ref  => catalog_accuracy($cat->{catalog}),
                               obsid    => $obsid,
                               detector => uc $cat->{detector},
                               det_ref  => $det_ref{ $cat->{detector} },
                               dr => $cat->{dr},
                               map { $_ => $cat->{$_} } qw(tstart dy dz catalog version date_obs
                                                           c_id x_id r_angle near_neighbor_dist double_id
                                                           status_id process_status grating sim_z net_counts snr
                                                           fids target x_ra x_dec c_ra c_dec),
                             };
        }
    }

    return @plot_data;
}

##*****************************************************************************
sub plot_all_offsets {
##*****************************************************************************
    my $file = shift;
    my $d = Ska::HashTable->new($file)->row("pos_ref < $opt{pos_ref_lim}");

    my @det_names = qw(ACIS-S ACIS-I HRC-S HRC-I);

    my $device = "$PLOT_DIR/offsets.ps/cps";	# color postscript file
    $device = '/xs';

    my $t = pdl $d->col("tstart");
    my $yr = $t / 86400. / 365.25 + 1998.0; # /  (something wrong with emacs perl mode)
    my ($x0, $x1) = minmax($yr);
    $x0 -= 0.1;
    $x1 += 0.1;
    my ($y0, $y1) = (-3, 3);

    # Make a HashTable object for each of the four detectors
    my @d_for_det = map { scalar $d->row("det_ref eq $_") } (0..3);

    my $plot_info = {data => \@d_for_det,
                     det_num => 0,
                     interactive => 1,
                     plot_type   => 'strip',
                     x0 => $x0,
                     x1 => $x1,
                     y0 => $y0,
                     y1 => $y1,
                     x0_lim => $x0,
                     x1_lim => $x1,
                     y0_lim => $y0,
                     y1_lim => $y1,
                     charsize => 1.4,
                     symbol_size => 1.0,
                     xsize => 12,
                     ysize => 8,
                     device => '/xs',
                     det_names => \@det_names,
                     plot_dr  => 0,
                     plot_bad => 0,
                     title => $title{dy_dz},
                    };

    if ($opt{batch}) {
        foreach my $det (0..3) {
            # Make strip plot of Y,Z offsets vs. time
            my $strip_plot_info = { %$plot_info,
                                    device      => "$PLOT_DIR/offsets-$det_names[$det].ps/cps",
                                    plot_type   => 'strip',
                                    det_num     => $det,
                                    interactive => 0,
                                    charsize    => 2.0,
                                    symbol_size => 1.5,
                                  };
            plot_offsets( $strip_plot_info );

            # Make histogram of radial offsets
            my $hist_plot_info = { %$plot_info,
                                   device      => "$PLOT_DIR/offsets-$det_names[$det]-hist.ps/cps",
                                   det_num     => $det,
                                   interactive => 0,
                                   charsize    => 2.0,
                                   symbol_size => 2.0,
                                   plot_type   => 'r_histogram',
                                   x0 => 0,
                                   x1 => 1.2,
                                   y0 => 0,
                                   y1 => 1,
                                   x0_lim => 0,
                                   x1_lim => 1.2,
                                   y0_lim => 0,
                                   y1_lim => 1,
                                   xsize => 8,
                                   ysize => 6,
                                   plot_dr => 1,
                                 };
            plot_offsets( $hist_plot_info );

            # Make histogram of radial offsets
            $hist_plot_info = { %$plot_info,
                                   device      => "$PLOT_DIR/offsets-$det_names[$det]-xyhist.ps/cps",
                                   det_num     => $det,
                                   interactive => 0,
                                   charsize    => 2.0,
                                   symbol_size => 2.0,
                                   plot_type   => 'xy_histogram',
                                   x0 => -1.5,
                                   x1 => 1.5,
                                   y0 => 0,
                                   y1 => 1,
                                   x0_lim => -1.5,
                                   x1_lim => 1.5,
                                   y0_lim => 0,
                                   y1_lim => 1,
                                   xsize => 8,
                                   ysize => 6,
                                 };
            plot_offsets( $hist_plot_info );

            # Make cumulative histogram of x,y offsets
            $hist_plot_info = { %$plot_info,
                                   device      => "$PLOT_DIR/offsets-$det_names[$det]-cxyhist.ps/cps",
                                   det_num     => $det,
                                   interactive => 0,
                                   charsize    => 2.0,
                                   symbol_size => 2.0,
                                   plot_type   => 'cxy_histogram',
                                   x0 => -1,
                                   x1 => 1,
                                   y0 => 0,
                                   y1 => 1,
                                   x0_lim => -1,
                                   x1_lim => 1,
                                   y0_lim => 0,
                                   y1_lim => 1,
                                   xsize => 8,
                                   ysize => 6,
                                 };
            plot_offsets( $hist_plot_info );
        }
    } else {
        while (plot_offsets($plot_info)) {};
    }
}

##*****************************************************************************
sub plot_offsets {
##*****************************************************************************
    my $pi = shift;		# Plot Info

    unless (defined $pi->{win}) {
        $pi->{win} = PDL::Graphics::PGPLOT::Window->new(Device => $pi->{device},
                                                        WindowXSize => $pi->{xsize},
                                                        WindowYSize => $pi->{ysize},
                                                       );
        $pi->{win}->ctab(lut_data('rainbow2')); # set up color palette to 'idl5'
    }

    my $win = $pi->{win};

    # Draw the new plot window
    $win->env($pi->{x0}, $pi->{x1}, $pi->{y0}, $pi->{y1},
              {charsize => $pi->{charsize}}
             );

    # Define the data.  Each var is an array ref
    my @points = get_plot_points($pi);

    if ($pi->{plot_type} eq 'r_histogram') {
        my ($p) = collate(@points);
        my $n_points = @{$p->{x}};
        my ($x, $y) = hist pdl($p->{y});
        my $cy = cumusumover $y;
        $cy /= $cy($cy->nelem-1);
        $win->bin($x, $cy);
        $win->label_axes('Radial offset (arcsec)',
                         'Cumulative fraction',
                         $pi->{det_names}[$pi->{det_num}] . " ($n_points points)",
                         {charsize=>$pi->{charsize}});
        $win->hold();
        plot_hist_lim($win, $x, $cy, 0.9, '90% limit');
        plot_hist_lim($win, $x, $cy, 0.68, '68% limit');
    } elsif ($pi->{plot_type} eq 'vector') {
        $win->label_axes('Y angle (arcmin)', 'Z angle (arcmin)', $pi->{det_names}[$pi->{det_num}],
                         {charsize=>$pi->{charsize}});
        $win->hold();
        foreach my $p (@points) {
        }
    } elsif ($pi->{plot_type} eq 'strip') {
    # Collate (by color and symbol) the points and plot
        for my $p (collate(@points)) {
            $win->points($p->{x}, $p->{y}, {color => $p->{color}, symbol => $p->{symbol}, charsize=>$pi->{symbol_size}});
        }

        $win->label_axes("Year",
                         "Offset (arcsec)",
                         $pi->{det_names}[$pi->{det_num}] . " : $pi->{title}",
                         {charsize=>$pi->{charsize}});

        # Draw reference dashed lines to show -2..+2 arcsec offsets
        foreach (-2,-1,0,1,2) {
            $win->line([$pi->{x0},$pi->{x1}], [$_, $_], {linestyle => 'dotted'});
        }
    } elsif ($pi->{plot_type} eq 'cxy_histogram') {
        # Collate (by color and symbol) the points and plot
        my %color_map = ('blue' => 'dz',
                         'red' => 'dy');
        for my $p (collate(@points)) {
            my ($x, $y) = hist(pdl($p->{'y'}), -3, 3, 0.02);
            my $cy = cumusumover $y;
            $cy /= $cy($cy->nelem-1);
            $win->bin($x, $cy, {color => $p->{color}});
            printf("%-6s : %-3s : %6.2f\n",
                   $pi->{det_names}->[$pi->{det_num}],
                   $color_map{$p->{color}},
                   median(pdl($p->{'y'}))
                  );
        }

        $win->label_axes("Offset (arcsec)",
                         "Distribution",
                         $pi->{det_names}[$pi->{det_num}] . " : $pi->{title}",
                         {charsize=>$pi->{charsize}});
        # Draw reference dashed lines to show 50%
        $win->line([$pi->{x0},$pi->{x1}], [0.5, 0.5], {linestyle => 'dotted'});

    } elsif ($pi->{plot_type} eq 'xy_histogram') {
        # Collate (by color and symbol) the points and plot
        for my $p (collate(@points)) {
            my ($x, $y) = hist(pdl($p->{'y'}), -3, 3, 0.1);
            $y  /= max($y);
            $win->bin($x, $y,  {color => $p->{color}});
        }

        $win->label_axes("Offset (arcsec)",
                         "Distribution",
                         $pi->{det_names}[$pi->{det_num}] . " : $pi->{title}",
                         {charsize=>$pi->{charsize}});
    }

    # If batch mode or histogram or vector, just close window and return
    if ($pi->{plot_type} =~ /histogram|vector/ or not $pi->{interactive}) {
        $win->close;
        delete $pi->{win};
        return;
    }

    # Get key/mouse inputs and perform actions
    my $ch = '';
    ($pi->{x}, $pi->{y}, $ch) = $win->cursor({# Type => 'CrossHair',
                                              Xref => $pi->{xref},
                                              Yref => $pi->{yref},
                                             });
    if ($ch eq 'A') {
        zoom_in($pi, $win);
    } elsif ($ch eq 'p' or $ch eq 'X') {
        pan_out($pi);
    } elsif ($ch =~ /\A [1234] \Z/x) {
        $pi->{det_num} = $ch-1;
    } elsif ($ch =~ /[Dics ]/) {
        inspect_data_point($pi, $ch, @points);
    } elsif ($ch =~ /b/) {
        $pi->{plot_bad} = not $pi->{plot_bad};
    } elsif ($ch =~ /r/) {
        $pi->{plot_dr} = not $pi->{plot_dr};
        $pi->{title} = $pi->{plot_dr} ? $title{dr} : $title{dy_dz};
        $pi->{y0} = $pi->{plot_dr} ? 0 : -3;
        $pi->{y0_lim} = $pi->{plot_dr} ? 0 : -3;
    } elsif ($ch eq 'h') {
        print <<COMMANDS;
Commands:
  Left-click   : Zoom (left-click again to define other corner)
  Middle-click : Info about point at cursor
   or <space>
  Right-click  : Pan
   or 'p'
  1 2 3 4      : Change detector (ACIS-S, ACIS-I, HRC-S, HRC-I)
  i            : Make image of point at cursor (gets evt data from archive)
  c            : Clean files associated with point at cursor
  s            : Change status_id for point (set as extended etc)
  b            : Toggle whether to plot 'bad' points (with status_id != 0)
  r            : Toggle plotting dy/dz or dr
  h            : List this help
  q            : Quit
COMMANDS
;
    }

    $pi->{xref} = $pi->{x};
    $pi->{yref} = $pi->{y};


    return ($ch ne 'q');
}

##*****************************************************************************
sub plot_hist_lim {
##*****************************************************************************
    my $win = shift;
    my $x = shift;
    my $cy = shift;
    my $lim = shift;
    my $txt = shift;
    my $offset = shift || -0.01;
    my $ok = which ($cy > $lim);
    my $x0 = $x($ok(0));
    my $y0 = $cy($ok(0));
    my $x20 = $x0->append($x0);
    $win->line($x20, pdl(-1,$y0->at(0)), {linestyle=>'DASHED'});
    $win->line(pdl(-1,$x0->at(0)), pdl($lim,$lim), {linestyle=>'DASHED'});
    $win->text($txt, $x0->at(0) + $offset, 0.1, {Angle=>90});
}

##*****************************************************************************
sub collate {
# Take a single array of points with different symbols and colors and collect into
# arrays for each symbol and color
##*****************************************************************************
    my %symbol = ();
    my %color  = ();
    my %coll = ();
    my @points = @_;
    my @out;

    # Find all the symbols and colors
    for my $p (@points) {
        $symbol{$p->{symbol}} = 1;
        $color{$p->{color}} = 1;
    }

    # Do the collection.  Not the most efficient algorithm (repeated greps) but
    # performance is not an issue
    for my $s (keys %symbol) {
        for my $c (keys %color) {
            my @p_ok = grep { $_->{color} eq $c and $_->{symbol} eq $s } @points;
            next unless @p_ok;
            my %out = ();
            for my $f (qw(x y data)) {
                $out{$f} = [ map { $_->{$f} } @p_ok ];
            }
            $out{symbol} = $s;
            $out{color} = $c;
            push @out, \%out;
        }
    }

    return @out;
}

##*****************************************************************************
sub get_plot_points {
##*****************************************************************************
    my $pi = shift;
    my @point;

    my %symbol = ( dy => 13,
                   dz => 16,
                   dr => 16,
                 );
    my %color = ( 'ok'  => { dy => 'red',
                             dz => 'blue',
                             dr => 'green',
                           },
                  'bad' => { dy => 'yellow',
                             dz => 'green',
                             dr => 'blue',
                           },
                );
    my %off_axis = ( dy => 'y_angle',
                     dz => 'z_angle',
                     dr => 'r_angle'
                   );

    # Define axes and status vals to plot
    my @plot_axis   = $pi->{plot_dr}  ? qw(dr)     : qw(dy dz);
    my @plot_status = $pi->{plot_bad} ? qw(ok bad) : qw(ok);

    my $d_det = $pi->{data}->[$pi->{det_num}];
    $d_det->reset();
    while (my %d = $d_det->row) {
        for my $status (@plot_status) {
            next if ($status eq 'bad' and $d{status_id}==0);
            next if ($status eq 'ok'  and $d{status_id}!=0);
            for my $axis (@plot_axis) {
                push @point, { x      => $d{tstart} / 86400. / 365.25 + 1998.0,
                               y      => $d{$axis},
                               oaa    => $d{ $off_axis{$axis} },
                               color  =>  $color{$status}{$axis},
                               symbol => $symbol{$axis},
                               data   => \%d,
                             };
            }
        }
    }

    return @point;
}

##*****************************************************************************
sub inspect_data_point {
##*****************************************************************************
    my ($pi, $ch, @points) = @_;
    my $src = find_nearest($pi, @points);

    printf("Record nearest %.5f, %.2f: Obsid=%d x_id=%d c_id=%d\n",
           $pi->{x}, $pi->{y}, $src->{obsid}, $src->{x_id}, $src->{c_id});
    $Data::Dumper::Sortkeys = 1;
    print Dumper $src;

    # Create and display image of source
    if ($ch eq 'i') {
        my $img_jpeg = make_source_image($src);
        system("display $img_jpeg &");
    }
    # Clean files associated with the data point
    elsif ($ch eq 'c') {
        print "Clean: (e)vent file, (i)mage files, (a)ll :";
        (undef, undef, $ch) = $pi->{win}->cursor();
        print "$ch\n";
        my $dir = io(obs_dir($src->{obsid}));
        $dir->chdir;
        my $fileglob;
        $fileglob = "acis*evt2.fits* hrc*evt2.fits*" if $ch eq 'e';
        $fileglob = "img_*"                          if $ch eq 'i';
        $fileglob = "*"                              if $ch eq 'a';
        my @files = glob($fileglob);
        if (@files) {
            print "Removing: ", join("\n",@files), "\n";
            print "OK (N/y)? ";
            (undef, undef, $ch) = $pi->{win}->cursor();
            print "$ch\n";
            unlink @files if ($ch eq 'y');
        } else {
            print "No matching files\n";
        }
    }
    # Set status_id
    elsif ($ch eq 's') {
        print "Set status_id for x_id=$src->{x_id}.  Current status_id=$src->{status_id}\n";
        my @status = sql_fetchall_array_of_hashref($dbh,
                                                   "select * from astromon_status order by id");
        foreach (@status) {
            printf("%-3d: %s\n", $_->{id}, $_->{status});
        }
        print "Value: ";
        (undef, undef, $ch) = $pi->{win}->cursor();
        print "$ch\n";

        if (grep { $ch eq $_->{id} } @status) {
            print "Setting status_id to $ch\n";
            my $update = <<END_UPDATE;
UPDATE astromon_xray_src SET status_id=$ch
 WHERE obsid=$src->{obsid}
   AND id=$src->{x_id}
END_UPDATE
            sql_do($dbh, $update) or die "Update failed: ", $dbh->errstr;
        }
    }


}

##*****************************************************************************
sub find_nearest {
##*****************************************************************************
    my $pi = shift;
    my @points = @_;

    my $x = pdl( map { $_->{x} } @points );
    my $y = pdl( map { $_->{y} } @points );

    my $delx = ($x - $pi->{x}) / ($pi->{x1}-$pi->{x0});
    my $dely = ($y - $pi->{y}) / ($pi->{y1}-$pi->{y0});
    my (undef, undef, $min_index) = minmaximum( $delx**2 + $dely**2 );

    return $points[ sclr($min_index) ]->{data};
}

##*****************************************************************************
sub zoom_in {
##*****************************************************************************
    my $pi = shift;
    my $win = shift;
    my $ch;

    my $x0 = $pi->{x};
    my $y0 = $pi->{y};

    ($pi->{x}, $pi->{y}, $ch) = $win->cursor({Type => 'Rectangle',
                                              Xref => $x0,
                                              Yref => $y0,
                                             });
    return if ($ch ne 'A');

    $pi->{x0} = $x0;
    $pi->{y0} = $y0;

    for my $a (qw(x y)) {
        $pi->{"${a}1"} = $pi->{${a}};
        if ($pi->{"${a}0"} > $pi->{"${a}1"}) {
            my $tmp = $pi->{"${a}0"};
            $pi->{"${a}0"} = $pi->{"${a}1"};
            $pi->{"${a}1"} = $tmp;
        }
        if ($pi->{"${a}0"} == $pi->{"${a}1"}) {
            $pi->{"${a}0"} = $pi->{"${a}0_lim"};
            $pi->{"${a}1"} = $pi->{"${a}1_lim"};
        }
        $pi->{$a} = ($pi->{"${a}0"} + $pi->{"${a}1"})/2;
    }
}


##*****************************************************************************
sub pan_out {
##*****************************************************************************
    my $pi = shift;

    for my $a (qw(x y)) {
        my $mid = ($pi->{"${a}0"} + $pi->{"${a}1"})/2;
        my $dist = ($pi->{"${a}1"} - $pi->{"${a}0"});
        $pi->{"${a}0"} = $mid - $dist;
        $pi->{"${a}1"} = $mid + $dist;
        $pi->{"${a}0"} = $pi->{"${a}0_lim"} if ($pi->{"${a}0"} < $pi->{"${a}0_lim"});
        $pi->{"${a}1"} = $pi->{"${a}1_lim"} if ($pi->{"${a}1"} > $pi->{"${a}1_lim"});
    }
}

##*****************************************************************************
sub unique_sources {
# Get the "best" corresponding catalog source for each X-ray source
##*****************************************************************************
    my @db = @_;
    my @src = ();
    my $match;

    # Get all the unique X-ray IDs
    my @x_ids = List::MoreUtils::uniq map { $_->{x_id} } @db;

    # For each unique X-ray ID, find all the corresponding xcorr db entries
    # and find the "best" catalog among those
    foreach my $x_id (@x_ids) {
        my @db_for_x_id = grep { $_->{x_id} == $x_id } @db;
        my $src;
        foreach my $db (@db_for_x_id) {
            $src = $db if (not defined $src
                           or $opt{all}
                           or $db->{version} > $src->{version}
                           or catalog_accuracy($db->{catalog}) < catalog_accuracy($src->{catalog})
                           or $db->{dr} < $src->{dr}
                          );
        }
        push @src, $src;
    }

    return @src;
}

##*****************************************************************************
sub catalog_accuracy {
# Return a numerical code ranking accuracy
#  0-9   :   < 0.1 arcsec
#  10-19 :   < 0.5 arcsec
#  20    :   > 0.5 arcsec
##*****************************************************************************
    my $catalog = shift;
    my %accuracy = ( ICRS => 0,
                     Tycho => 1,
                     SIMBAD_high => 2,
                     CELMON => 3,
                     ASTROMON => 4,
                     SDSS => 10,
                     '2MASS' => 12,
                     SIMBAD_med => 14,
                     USNO => 20,
                   );

    my $accuracy = 20;
    for (keys %accuracy) {
        $accuracy = $accuracy{$_} if ($catalog =~ /$_/i);
    }

    return $accuracy;
}

##****************************************************************************
sub make_source_image {
# Use PGPLOT to make color image map of source and immediate surrounding
# and the source extraction / background regions
##****************************************************************************
    # Define various file names for use later
    my $src = shift;
    my $obsid = $src->{obsid};
    my $dir = io(obs_dir($obsid));
    $dir->chdir;
    (my $instrument = lc $src->{detector}) =~ s/-.+//;

    my $img_ps = sprintf("img_obs%d_x%d_c%d.ps", $src->{obsid}, $src->{x_id}, $src->{c_id});
    (my $img_fits = $img_ps) =~ s/\.ps\Z/.fits/;
    (my $img_jpeg = $img_ps) =~ s/\.ps\Z/.jpg/;

    return io($img_jpeg)->absolute->pathname if (-r $img_jpeg);

    print "Getting evt2 data for obsid $obsid: ";
    my ($evt2) = get_archive_files(obsid     => $obsid,
                                   prod      => $instrument . '2{evt2}',
                                   file_glob => $instrument . '*evt2.fits*',
                                   dir       => "$dir",
                                  );
    die "Could not get evt2 file\n" unless $evt2;
    print "$evt2\n";

    # Set global hash of event files to allow option of deleting on program exit
    $evt2_files{$evt2} = io($evt2)->absolute->pathname;

    # Weird problem with remote arc5gl causing dismount of /soft/ciao/CALDB.
    # See: http://jeeves.cfa.harvard.edu/CXCAspect/Aspect/RemoteArcFiveGlProblem
    while (not -d "$ENV{CALDB}/") {
        print "Waiting to see $ENV{CALDB}/\n";
        sleep 2;
    }

    print "Running dmcoords...\n";
    run_tool("punlearn dmcoords");
    confess "Dmcoords error" unless my $DMcoord = Chandra::Tools::dmcoords->new( $evt2 );

    my ($out) = $DMcoord->coords( cel => ($src->{x_ra}, $src->{x_dec}));
    my $x_x = $out->{sky}{'x'};
    my $x_y = $out->{sky}{'y'};

    ($out) = $DMcoord->coords( cel => ($src->{c_ra}, $src->{c_dec}));
    my $c_x = $out->{sky}{'x'};
    my $c_y = $out->{sky}{'y'};

    my $bin = $instrument eq 'acis' ? 0.5 : 2;
    my $skypix_per_arcsec = $bin / 0.25; # This is tied to definition of $bin
    my $sz2 = $bin * 50;
    my $energy = $instrument eq 'acis' ? '[energy=500:8000]' : '';
    run_tool(sprintf("dmcopy '%s[bin x=%d:%d:%f,y=%d:%d:%f]$energy' $img_fits clobber=yes",
                     $evt2, $x_x-$sz2, $x_x+$sz2, $bin, $x_y-$sz2, $x_y+$sz2, $bin),
             { loud=>1 }
            );

    my $win = PDL::Graphics::PGPLOT::Window->new(device => "$img_ps/vcps",
                                                 size => 5,
                                                 unit => 'inch', # Inches
                                                 hold => 1,
                                                 axis =>  'box',  # Make plain (empty) axes
                                                );
    $win->ctab(lut_data('heat', 1));
    my $img = rfits($img_fits);

    # Force fits_imag to use Physical (sky) WCS coords instead of (RA,Dec)
    # which is the default
    for my $axis (1, 2) {
        for my $hdrkey (qw(CTYPE CRVAL CRPIX CDELT)) {
            $img->hdr->{"${hdrkey}${axis}"} = $img->hdr->{"${hdrkey}${axis}P"};
        }
    }
    $win->fits_imag($img, { itf => 'log',
                            DrawWedge => 0,
                            xtitle => ' ',
                            ytitle=> ' ',
                          });
    $win->hold(1);
    my %figure_opt = ( linewidth => 4,
                       filltype => 'outline',
                       color => 'green',
                     );
    # make a circle with radius 1"
    $win->circle($x_x, $x_y, $skypix_per_arcsec, {%figure_opt});
    $win->line([$c_x-$skypix_per_arcsec/2, $c_x+$skypix_per_arcsec/2],
               [$c_y, $c_y],
               { %figure_opt });
    $win->line(	[$c_x, $c_x],
                [$c_y-$skypix_per_arcsec/2, $c_y+$skypix_per_arcsec/2],
                { %figure_opt },
              );

    $win->close();

    # Convert to jpeg by using the eps2png perl script
    run_tool("ps2any $img_ps -size 400",
             { loud=>1 });
    unlink $img_ps if (-e $img_ps);

    return io($img_jpeg)->absolute->pathname;
}

sub obs_dir {
    my $obsid = shift;
    my $dir = sprintf("$ASTROMON_DATA/obs%02d/%d", $obsid/1000, $obsid);
    io($dir)->mkpath unless io($dir)->exists;
    return $dir;
}

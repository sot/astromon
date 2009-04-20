#!/usr/bin/env /proj/sot/ska/bin/perl

use warnings;
use strict;

##########################################################################################
#
# Monitor celestial location
#
# Get obsids (could be spec'd on cmd line, or via standard status tables with date range)
# Remove observations of SNR and solar system objects
# Read existing cel_loc table
#
# For each obsid
#   Extract events file, source file if possible
#   if (grating)
#     Find 0-order source with celldetect and create src2 file
#   Set sources to values in src2 file (check for point-like-ness)
#   Remove piled sources or SNR < SNR_MIN
#   for each (@sources)
#     cross correlate with catalogs and find best match
#     if (matched) 
#       Update cel_loc table
#       Make jpg of source (log scale) (dmimg2jpg)
#   Update cel_loc table
#
# Update cel_loc plots
##########################################################################################

use Getopt::Long;
use App::Env qw( CIAO );

use Ska::Convert qw(time2date);
use Ska::Process qw(get_params get_archive_files);
use Ska::CatQuery qw(get_seqnum_info);
use Ska::Message;
use IO::All;
use Data::Dumper;
use DBI;
use DBD::Sybase;
use POSIX qw(strftime);

our $ASTROMON_SHARE = "$ENV{SKA}/share/astromon";
our $ASTROMON_DATA  = "$ENV{SKA}/data/astromon";
# Set up some constants

$ENV{UPARM} = $ASTROMON_DATA;
$| = 1;

##****************************************************************************
# MAIN PROGRAM
##****************************************************************************

our %par = get_options();	# Parse command line options

my $log_dir = io("$ASTROMON_DATA/logs/get_cat_obs_data")->mkpath;
our $log = Ska::Message->new(file => "$log_dir/" . strftime("%Y-%m-%d_%H:%M", localtime));

my $date = localtime;
my $banner = "----------- Starting get_cat_obs_data.pl at $date ------------";
$log->message("");
$log->message('-'x length $banner, "\n$banner\n", '-'x length $banner);
$log->message("");
print_options();

my @obsid = get_obsids();	# Extract the list of obsids

our %Category_ID_Map = (
		     'Normal Stars and WD'      => 10,
		     'WD Binaries and CVs'      => 20,
		     'BH and NS Binaries'       => 30,
		     'Normal Galaxies'          => 40,
		     'Active Galaxies and Quasars' => 50,
		     'Extragalactic Diffuse Emission & Surveys' => 60,
		     'Galactic Diffuse Emission & Surveys' => 70,
		     'Solar System and Misc'    => 100,
		     'SN, SNR, and Isolated NS' => 110,
		     'Clusters of Galaxies'     => 120,
		     'NO MATCH'                 => 200,
		    );
$Category_ID_Map{lc($_)} = $Category_ID_Map{$_} for keys %Category_ID_Map;

# Get sybase database handle for aca database for user=aca_ops
my $dbh = get_dbh();

my %astromon_table_def = eval scalar io("$ASTROMON_DATA/astromon_table_defs")->slurp;

OBSID: foreach my $obsid (@obsid) {
    $date = localtime;
    $log->message("********* PROCESSING OBSID $obsid at $date ***********");

    my $obs = Obs->new(obsid => $obsid);

    my %table = (astromon_obs      => {},
		 astromon_xray_src => [],
		 astromon_cat_src  => [],
		);			# Final tables that will get put to database

    # Check for an existing entry in database for this observation
    my $astromon_obs = get_astromon_obs_from_db($dbh, $obsid);

    # Don't redo cross correlation, don't log message
    if (defined $astromon_obs->{process_status} and not $par{force}) {	
	$log->message("  Already processed");
	next OBSID;
    }

    # These routines throw exceptions for any problems or stop signals in processing
    eval {
	my %td;			# Table definition hashes
	# Set a working directory, make it, and cd there
	$obs->work_dir(sprintf("$par{dirs}/obs%02d/%d", $obsid/1000, $obsid));
	io($obs->work_dir)->mkpath;
	chdir $obs->work_dir;
    
	# Clean files if needed
	$obs->preclean($par{preclean_files});

	# Get the category (e.g. Extragalactic) for obs and map to numerical value for database
	# (Should be done self-consistently with DB table, foreign keys etc!)
	my %info = get_seqnum_info($obs->obspar->{seq_num});
	$obs->category_id( $Category_ID_Map{lc($info{categ})} || $Category_ID_Map{'NO MATCH'} );

	# Set all the fields in the astromon_obs table
	%td = grep {not ref} @{$astromon_table_def{astromon_obs}}; # How to coerce table to hash in one step??
	foreach (sort keys %td) {$table{astromon_obs}->{$_} = $obs->$_;} 

	# Get the X-ray source list for this obsid
	%td = grep {not ref} @{$astromon_table_def{astromon_xray_src}}; 
	my $src2_data = Data::ParseTable::parse_table($obs->src2); 

	my @xray_sources = map { XraySource->new(obs  => $obs, src2 => $_) } @$src2_data;
	foreach my $xray_source (@xray_sources) {
	    # Tell source about the other sources
	    $xray_source->xray_sources([ grep { $xray_source != $_ } @xray_sources ]); 
	    push @{$table{astromon_xray_src}}, { map { $_ => $xray_source->$_ } sort keys %td };
	}
	    
	# Get counterpart information from various catalogs
	$table{astromon_cat_src} = $obs->get_catalog_source_list();

	$table{astromon_obs}->{process_status} = 'OK';
    };

    if ($@) {
	# Observation not fully processed for some reason.  Could be benign (wrong
	# obs type) or serious (couldn't get a needed file)
	chomp $@;
	my $error = ($@ =~ s/^OK: //) ? '' : 'ERROR: ';
	$log->message("  ${error}$@\n");

	# Explicitly set obsid and process status for this failure case
	$table{astromon_obs}->{obsid} = $obsid;
	$table{astromon_obs}->{process_status} = $@;
    }
	
    # Finally store the available data, writing over any existing entries.
    for my $table_name (qw(astromon_obs astromon_xray_src astromon_cat_src)) {
	overwrite_table_rows($dbh,         # database handle
			     $table_name,  # database table name
			     {'obsid' => $obsid},  # key columns and value for dropping existing records
			     $astromon_table_def{$table_name},  # array ref of table definition
			     $table{$table_name}  # actual data (either array ref or hash ref)
			    );
    }

    $obs->clean_up($par{clean});
    $obs->gzip_fits();
}

$dbh->disconnect();

##****************************************************************************
sub overwrite_table_rows {
##****************************************************************************
    my $dbh = shift;       # database handle
    my $table_name = shift;	# database table name
    my $delete_keys = shift;	# key columns and values for dropping existing records
    my $table_def  = shift;	# array ref of table definition
    my $table_data = shift;	# Actual data (either array ref or hash ref)
				#  that may contain extra key/val pairs

    $log->message("Writing table $table_name");

    # Force table data into an array of table rows (each one being a hashref)
    my @table_data = ref($table_data) eq 'ARRAY' ? @$table_data : ($table_data);

    # Get the columns from table definition (which is an array ref of col=>def pairs)
    my @table_cols = keys %{ {grep {not ref} @$table_def} };

    my $where = join(' AND ', map { "$_=$delete_keys->{$_}" } keys %$delete_keys);
    my $do = "DELETE FROM $table_name WHERE $where";
    $log->message($do);
    sql_do($dbh, $do);

    foreach my $data (@table_data) {
	# Insert new row of data, passing only the cols that are in the DB table
	insert_hash($dbh, $table_name, { map { $_ => $data->{$_} } @table_cols })
	  or die "Failed to insert a row into $table_name";
    }
}


##****************************************************************************
sub get_options {
##****************************************************************************
    my %par = get_params("astromon.par",
			 "tstart=s",
			 "tstop=s",
			 'caldb=s',
			 "month!",
			 'preclean_files=s',
			 'force',
			 'nodb',
			);

    die unless (%par);

    if ($par{month}) {
	my $time_now = $par{tstart} || time2date(time,1); # Current time
	my $date_now = `time_convert.pl $time_now fits`;
	my ($year_now,$month_now) = ($date_now =~ /(\d\d\d\d)-(\d\d)/);
	my $month_last = sprintf("%02d", ($month_now == 1) ? 12 : $month_now-1);
	my $year_last = ($month_now == 1) ? $year_now-1 : $year_now;
	$par{tstart} = "$year_last-$month_last-01T00:00:00";
	$par{tstop}  = "$year_now-$month_now-01T00:00:00";
    }

    $par{dirs} = "$ASTROMON_DATA/$par{dirs}" unless io($par{dirs})->is_absolute;

    return %par;
}

##****************************************************************************
sub get_dbh {
# Accessor returning (and possibly creating) global sybase handle
##****************************************************************************
    my $sybase_pwd < io("$ENV{SKA}/data/aspect_authorization/sybase-aca-aca_ops");
    chomp $sybase_pwd;

    ##-- Make database connection 2: sybase 
    my $server 	= "server=sybase";
    my $db     	= "database=aca";
    my $db_user = "aca_ops";
    my $dbh = DBI->connect("DBI:Sybase:$server;$db", $db_user, $sybase_pwd,
			    { PrintError => 0,
			      RaiseError => 1
			    }
			   );
    die "Couldn't connect: $DBI::errstr" unless ($dbh);

    return $dbh;
}

##****************************************************************************
sub print_options {
##****************************************************************************
    $log->message("COMMAND LINE PARAMETERS");
    foreach (sort keys %par) {
	$log->message(sprintf("  %-16s = %s", $_, $par{$_}));
    }
    $log->message("");
}

##****************************************************************************
sub get_obsids {
##****************************************************************************
    # if the argument is a file containing a list of obsids, read it, 
    # take the first column, and accept only numeric entries
    if (@ARGV && $ARGV[0] =~ /^@/) {
	$ARGV[0] =~ s/@//;
	@ARGV = grep /\d+/, `cat $ARGV[0]`;
	map {chomp; ($_) = split} @ARGV;
	@ARGV = grep /^\d+$/, @ARGV;
    }

    # The argument list can be a set of directories in the form obs<obsid>
    map {s/.*obs//; s/\///} @ARGV;

    # If tstart/tstop are defined, then go the archive to pull out all the
    # obsids in the specified time range
    if (exists $par{tstart}) {
	$par{tstart} = time2date(time+$par{tstart}*86400,1) if ($par{tstart} =~ /^-/);
	if ($par{tstop}) {
	    $par{tstop} = time2date(time+$par{tstop}*86400,1) if ($par{tstop} =~ /^-/);
	} else {
	    $par{tstop} = time2date(time,1);
	}

	$log->message("Processing obsids from $par{tstart} to $par{tstop}\n");

	my $obspar_dir = io("$ASTROMON_DATA/obspar");
	$obspar_dir->rmtree if $obspar_dir->exists;
	$obspar_dir->mkpath;

	my @obsfiles =  get_archive_files(tstart    => $par{tstart},
					  tstop     => $par{tstop},
					  prod      => "obspar",
					  file_glob => "axaf*obs0a.par*",
					  dir       => "$obspar_dir",
					  gunzip    => 0,
                                          version   => [''],
					  loud      => 0);
	undef $Ska::Process::arc5gl;
	$obspar_dir->rmtree if $obspar_dir->exists;

	push @ARGV, grep {$_ = $1+0 if (/axaff(\d+)/ && $1 < 50000)} @obsfiles;
    }

    $log->message("Processing following obsids: @ARGV");
    return @ARGV;
}

##****************************************************************************
sub get_astromon_obs_from_db {
##****************************************************************************
    my $dbh = shift;
    my $obsid = shift;
    my $astromon_obs;

    my @row = fetchall_array_of_hashref($dbh,
					"select * from astromon_obs where obsid=?", $obsid);
    if (@row == 0) {
	$astromon_obs = {};
    } elsif (@row == 1) {
	$astromon_obs = $row[0];
    } else {
	die "Database is corrupted: multiple entries in astromon_obs for obsid $obsid\n";
    }

    return $astromon_obs;
}

##************************************************************************---
  sub insert_hash {
##***************************************************************************
    my ($dbh, $table, $field_values) = @_;
    my @fields = sort keys %$field_values; # sort required
    my @values = @{$field_values}{@fields};
    my $sql = sprintf "insert into %s (%s) values (%s)",
        $table, join(",", @fields), join(",", ("?")x@fields);

    if ($par{nodb}) {
	print "Ignoring insert_hash: $sql\n";
	return 1;
    }

    my $sth = $dbh->prepare_cached($sql);
    return $sth->execute(@values);
}

#########################################################################
sub sql_do {
#########################################################################
    my $dbh = shift;
    my $sql = shift;
    if ($par{nodb}) {
	print "Ignoring dbh->do($sql)\n";
	return;
    }
    $dbh->do($sql) or die "$sql: " . $dbh->errstr;
}	


##****************************************************************************
sub fetchall_array_of_hashref {
##****************************************************************************
    my $dbh = shift;
    my $statement = shift;
    my @arg = @_;
    my @out;

    my $sth = $dbh->prepare($statement)
      or die "Bad SQL statement '$statement': " . $dbh->errstr;
    
    $sth->execute(@arg);
    while (my $row = $sth->fetchrow_hashref) {
	push @out, $row;
    }

    return @out;
 }

##************************************************************************---
sub calc_y_z_r {
##************************************************************************---
    my ($ra, $dec, $roll, $src_ra, $src_dec) = @_;
    my $r2d = 180/3.14159265358979;
    my $q = Quat->new($ra, $dec, $roll);
    my ($y_angle, $z_angle) = Quat::radec2yagzag($src_ra, $src_dec, $q);
    $y_angle *= $r2d * 3600;	# in arcsec
    $z_angle *= $r2d * 3600;	# in arcsec
    my $r_angle = sqrt( $y_angle**2 + $z_angle**2 );

    return ($y_angle, $z_angle, $r_angle);
}


##****************************************************************************
##****************************************************************************
package Obs;
##****************************************************************************
##****************************************************************************

use Ska::IO qw(read_param_file);
use Ska::Convert qw(dec2hms radec_dist hms2dec);
use Ska::CatQuery qw(get_simbad get_2mass);
use Ska::Process qw(get_archive_files);
use Data::ParseTable;
use RDB;
use IO::All;
use Data::Dumper;
use Ska::Web qw(get_url);
use XML::Simple;


use Class::MakeMethods::Standard::Hash ( 
					new   => 'new',
					scalar => [ qw( 
						       work_dir
						       obsid
						       process_status
						       category_id
						      )
						  ],
					array => 'cat',
					array => 'sources',
				       );

our $dbh;

##**************************************************************************--
sub obspar {
##****************************************************************************
    my $self = shift;

    return $self->{obspar} if defined $self->{obspar};

# Get obspar file, either from current directory or from archive

    $log->message("Getting obspar file");
    my ($obspar) = get_archive_files(obsid     => $self->obsid,
				     prod      => "obspar",
				     file_glob => "axaf*obs0a.par*",
				     dir       => $self->work_dir,
				  );
    die "Could not get an obspar file\n" unless defined $obspar;

    # Parse obspar
    my %obspar = read_param_file ($obspar);
    $obspar{ascdsver} =~ s/\A\s+|\s+\Z//g;

    die "Obsid mismatch ".$self->obsid." != $obspar{obs_id}\n" unless ($self->obsid == $obspar{obs_id});

    unless (exists $obspar{ra_targ} and exists $obspar{dec_targ}) {
	unless (exists $obspar{ra_nom} and exists $obspar{dec_nom}) {
	    die "No RA,DEC_NOM\n";
	}
	$obspar{ra_targ} = $obspar{ra_nom};
	$obspar{dec_targ} = $obspar{dec_nom};
    }

    if (exists $obspar{readmode} && $obspar{readmode} =~ /CONTINUOUS/i) {
	die "OK: Continuous clocking observation\n";
    }

    if ($obspar{tstop} - $obspar{tstart} > $par{max_time}*1000) {
	die sprintf("OK: Observation too long (%7.1f > %7.1f ksec)\n",
		    ($obspar{tstop} - $obspar{tstart})/1000., $par{max_time});
    }
    
    $log->message("Source from obspar file is $obspar{object}");

    return ($self->{obspar} = \%obspar);
}

##****************************************************************************
sub detector {
##****************************************************************************
    my $self = shift;
    ($self->{detector} = $self->obspar->{detector}) =~ tr/A-Z/a-z/
      unless defined $self->{detector};
    return $self->{detector};
}
    
##****************************************************************************
sub instrument {
##****************************************************************************
    my $self = shift;
    ($self->{instrument} = $self->detector) =~ s/-.+// unless defined $self->{instrument};
    return $self->{instrument};
}
    
##****************************************************************************
sub frame_time {
##****************************************************************************
    my $self = shift;
    unless (defined $self->{frame_time}) {
	if (($self->{frame_time} = ($self->obspar->{exptimea} || 1.0)) < 0.001) {
	    $self->{frame_time} = 1.0;
	}
    }
    return $self->{frame_time};
}

##****************************************************************************
sub extr_rad {
##****************************************************************************
    my $self = shift;
    $self->{extr_rad} = ($self->obspar->{grating} eq 'NONE' ? $par{extr_rad_image} : $par{extr_rad_grating})
      unless defined $self->{extr_rad};
    return $self->{extr_rad};
}

##****************************************************************************
sub ra_hms {
##****************************************************************************
    my $self = shift;
    unless (defined $self->{ra_hms}) {
	($self->{ra_hms}, $self->{dec_hms}) = dec2hms ($self->obspar->{ra_targ},
						       $self->obspar->{dec_targ});
    }
    return $self->{ra_hms};
}
	
##****************************************************************************
sub dec_hms {
##****************************************************************************
    my $self = shift;
    unless (defined $self->{dec_hms}) {
	($self->{ra_hms}, $self->{dec_hms}) = dec2hms ($self->obspar->{ra_targ},
						       $self->obspar->{dec_targ});
    }
    return $self->{dec_hms};
}
	
##****************************************************************************
sub target {
##****************************************************************************
    my $self = shift;
    $self->{target} = $self->obspar->{object} unless defined $self->{target};
    return $self->{target};
}
    
##****************************************************************************
sub grating {
##****************************************************************************
    my $self = shift;
    $self->{grating} = $self->obspar->{grating} unless defined $self->{grating};
    return $self->{grating};
}
    
##****************************************************************************
sub sim_z {
##****************************************************************************
    my $self = shift;
    $self->{sim_z} = sprintf("%.2f", $self->obspar->{sim_z}) unless defined $self->{sim_z};
    return $self->{sim_z};
}
    
##****************************************************************************
sub date_obs {
##****************************************************************************
    my $self = shift;
    $self->{date_obs} = $self->obspar->{'date-obs'} unless defined $self->{date_obs};
    return $self->{date_obs};
}
    
##****************************************************************************
sub tstart {
##****************************************************************************
    my $self = shift;
    $self->{tstart} = $self->obspar->{tstart} unless defined $self->{tstart};
    return $self->{tstart};
}
    
##****************************************************************************
sub ascdsver {
##****************************************************************************
    my $self = shift;
    $self->{ascdsver} = $self->obspar->{ascdsver} unless defined $self->{ascdsver};
    return $self->{ascdsver};
}

##****************************************************************************
sub ra {
##****************************************************************************
    my $self = shift;
    $self->{ra} = $self->obspar->{ra_nom} unless defined $self->{ra};
    return $self->{ra};
}

##****************************************************************************
sub dec {
##****************************************************************************
    my $self = shift;
    $self->{dec} = $self->obspar->{dec_nom} unless defined $self->{dec};
    return $self->{dec};
}

##****************************************************************************
sub roll {
##****************************************************************************
    my $self = shift;
    $self->{roll} = $self->obspar->{roll_nom} unless defined $self->{roll};
    return $self->{roll};
}

##****************************************************************************
sub cnt_rate {
##****************************************************************************
    my $self = shift;
    unless (defined $self->{cnt_rate}) {
	$self->{cnt_rate} = $self->counts / (($self->obspar->{tstop} - $self->obspar->{tstart})
					     / $self->frame_time);
    }
    return $self->{cnt_rate};
}

##****************************************************************************
sub version {
##****************************************************************************
    my $self = shift;
    return $self->{version} if defined $self->{version};

    my $version;
    if ($self->obspar->{ascdsver} =~ /R4CU(\d)UPD([\d\.]+)/) {
	$version = sprintf "%.3f", $1 + $2/100.0 ;
    } elsif ($self->obspar->{ascdsver} =~ /(\d+)\.(\d+)\.(\d+)/) {
	$version = sprintf "%.4f", $1 + $2/100.0 + $3/10000.0;
    } else {
	$version = "0.000";
    }
    return ($self->{version} = $version);
}

##****************************************************************************
sub get_source_list {
##****************************************************************************
    my $self = shift;
    
    # Read the cell detect source file, and then eliminate any sources

    my $src_data = Data::ParseTable::parse_table($self->src2); # Read src2 FITS file as array of hashes
    $self->sources( $src_data );
}

##****************************************************************************
sub src2 {
##****************************************************************************
    my $self = shift;

    # If src2 file not already defined, look in current directory
    if (not defined $self->{src2}) {
	$log->message("Getting src2 file");
	my ($src2) = get_archive_files(obsid     => $self->obsid,
				       prod      => $self->instrument . '2[*src2*]',
				       file_glob => "*src2.fits*",
				       gunzip    => 1,
				       dir       => $self->work_dir,
				      );

	# Not there, so go to archive to get evt2 file and then make src2
	if (not defined $src2) {
	    $log->message("Getting evt2 file");
	    my ($evt2) = get_archive_files(obsid     => $self->obsid,
					   prod      => $self->instrument . '2[*evt2*]',
					   file_glob => $self->instrument . '*evt2.fits*',
					   dir       => $self->work_dir,
					  );
	    die "Could not get evt2 file\n" unless $evt2;
	    $src2 = $self->make_src2_file($evt2);
	}
	$self->{src2} = $src2;
    }

    # Make sure we finally have a good src2 file
    die "Could not find or make a src2 file\n" unless defined $self->{src2};

    return $self->{src2};
}


##****************************************************************************
sub fidpr {
##****************************************************************************
    my $self = shift;

    if (not defined $self->{fidpr}) {
	$log->message("Getting fidpr file(s)");
	my @fidpr = get_archive_files(obsid     => $self->obsid,
				      prod      => "asp1[*fidpr*]",
				      file_glob => "pcad*fidpr*fits*",
				      version   => [4,3,2,1],
				      dir       => $self->work_dir,
				     );
	$self->{fidpr} = \@fidpr;
    }

    return $self->{fidpr};
}

##****************************************************************************
sub fids {
##****************************************************************************
    my $self = shift;
    my $COMMENT     = '#';

    return $self->{fids} if defined $self->{fids};

# Make a list of all fids used in observation.  There may be multiple
# intervals with different sets.   The routine would return something
# like '2,4,5:2,4' in this case.

    my %fids = ();		# All the fids used in observation 
    my @fidpr = @{$self->fidpr};
    $log->message("Reading fid props files @fidpr");
    foreach my $fidpr (@fidpr) {
	my $cmd = "dmlist \"$fidpr\[col id_num,id_status\]\" data,clean";
	my @lines = `$cmd`;
	my %fid = ();		# Fids used in one aspect interval
	foreach (@lines) {
	    next if (/^$COMMENT/);
	    my @a = split;
	    $fid{$a[0]} = 1 if ($a[1] eq 'GOOD');
	}
	$fids{join(',', sort {$a <=> $b} (keys %fid))} = 1; 
    }

    return (join ':', sort(keys %fids));
}

##****************************************************************************
sub make_src2_file {
##****************************************************************************
    my $self = shift;
    my $evt2 = shift;
    my $src_evt2 = "source_evt2.fits";
    my $src_src2 = "source_src2.fits";
    my $src_img  = "source_img.fits";
    my $img_bin  = ($evt2 =~ /\A hrc/x) ? "bin x=::8,y=::8" : "bin x=::1,y=::1";
    my $dmcopy;

# Extract a small radius around each celestial location source in field

    $log->message("Making src2 file");
    $log->message(sprintf("Extracting photons at %.5f,%.5f", $self->obspar->{ra_targ}, $self->obspar->{dec_targ}));
    $dmcopy = sprintf(qq{dmcopy "%s[events][(x,y)=circle(%s,%s,%f')]" %s clobber=yes},
		      $evt2, $self->ra_hms, $self->dec_hms, $self->extr_rad, $src_evt2);
    $log->message("$dmcopy");
    system($dmcopy);

    $dmcopy = "dmcopy \"${src_evt2}\[$img_bin\]\" $src_img clobber=yes";
    $log->message("$dmcopy");
    system($dmcopy);

# Find sources in the small field

    my $celldet = "celldetect $src_img $src_src2 thresh=$par{snr} clobber=yes";
    $log->message("$celldet");
    system($celldet);

    die "Failed to make src2 file\n" unless (-e $src_src2);

    return $src_src2;
}

##*************************************************************************---
sub get_catalog_source_list {
##****************************************************************************
    my $self = shift;

    $ENV{TY2_PATH} = $par{ty2_path};
    $ENV{UB1_PATH} = $par{ub1_path};

    my $cat;
    my @cat = ();
    my %data;
    my %obspar = %{$self->obspar};
    my $extr_rad = $self->extr_rad;

    my $t2000 = 63072000.0;	# Chandra time at 2000-01-01T00:00:00
    my $n_years_2000 = ($obspar{tstart} - $t2000) / (86400*365.25);
    my $year = sprintf "%.2f", 2000.0 + $n_years_2000;
    my $d2r = 3.14159265358979 / 180.;
    my $d2a = 3600.;

    # Read astromon_sources catalog

    $log->message("Reading ASTROMON catalog");
    my $rdb = new RDB "$ASTROMON_DATA/ASTROMON_table.rdb" or warn "Couldn't open $ASTROMON_DATA/ASTROMON_table.rdb\n";

    while($rdb && $rdb->read( \%data )) {
	my ($ra_hms, $dec_hms) = dec2hms($data{RA}, $data{Dec});
	next unless (radec_dist($data{RA},$data{Dec},$obspar{ra_targ},$obspar{dec_targ}) < $extr_rad/60.0);
	push @cat, {name    => $data{Name},
		    ra_hms  => $ra_hms,
		    dec_hms => $dec_hms,
		    ra      => $data{RA},
		    dec     => $data{Dec},
		    mag     => undef,
		    catalog => 'ASTROMON' };
    }
    
    # Read ICRS catalog

    $log->message("Reading ICRS catalog");
    my @cat_lines = `cat $ASTROMON_DATA/ICRS_tables`;
    foreach (@cat_lines) {
	next unless (/^J\d\d\d\d\d\d/);
	my @a = split;
	my ($ra, $dec) = hms2dec(@a[(2..7)]);
	next unless (radec_dist($ra, $dec, $obspar{ra_targ}, $obspar{dec_targ}) < $extr_rad/60.0);
	push @cat, {name    => $a[1],
		    ra_hms  => join(':',@a[(2..4)]),
		    dec_hms => join(':',@a[(5..7)]),
		    ra      => $ra,
		    dec     => $dec,
		    mag     => undef,
		    catalog => 'ICRS' };
    }

    my %cat_short_name = ('Tycho2'    => 'ty2',
			  'USNO-B1.0' => 'ub1');
    
    # Get 2MASS catalog sources

    if ($par{twomass}) {
	$log->message("Getting 2MASS catalog");
	my @cat_2mass = get_2mass($obspar{ra_targ}, $obspar{dec_targ}, $extr_rad*60);

	foreach my $data (@cat_2mass) {
	    my ($ra_hms, $dec_hms) = dec2hms($data->{ra}, $data->{dec});
	    push @cat, {name    => $data->{designation},
			ra_hms  => $ra_hms,
			dec_hms => $dec_hms,
			ra      => $data->{ra},
			dec     => $data->{dec},
			mag     => $data->{k_m}, # V mag
			catalog => '2MASS' };
	}
    }
    
    # Get Tycho and USNO catalogs

    $log->message("Reading Tycho2 and USNO B1 catalogs");
    for my $cat_name (keys %cat_short_name) {
	my $radius   = $extr_rad*60;
	my $scat_cmd = sprintf(qq{scat -n 1000 -y $year -c $cat_short_name{$cat_name} -r %.2f %s %s J2000},
			       $extr_rad*60, $self->ra_hms, $self->dec_hms);
	$log->message("$scat_cmd");
	my @cat_lines = `$scat_cmd`;

	foreach (@cat_lines) {
	    my @a = split;
	    my ($ra, $dec) = hms2dec($a[1], $a[2]);
	    next if ($a[4] > $par{ub1_mag_lim} and $cat_name eq 'USNO-B1.0');
	    push @cat, {name    => $a[0],
			ra_hms  => $a[1],
			dec_hms => $a[2],
			ra      => $ra,
			dec     => $dec,
			mag     => $a[4], # V mag
			catalog => $cat_name };
	}
    }

    # Get a coordinate from SIMBAD

    $log->message("Getting SIMBAD coordinate");
    my ($ra, $dec, $mag, $pm_ra, $pm_dec, $precision) = get_simbad($obspar{object});

    if (defined $ra and defined $dec) {
	# Apply proper motion correction since 2000.0 epoch
	if (defined $pm_ra and defined $pm_dec) {
	    $ra  += $pm_ra / $d2a / cos($dec*$d2r) * $n_years_2000;
	    $dec += $pm_dec / $d2a * $n_years_2000;
	}
	my ($rah, $decd) = dec2hms($ra, $dec);
	push @cat, {name    => $obspar{object},
		    ra_hms  => $rah,
		    dec_hms => $decd,
		    ra      => $ra,
		    dec     => $dec,
		    mag     => $mag || undef,
		    catalog => "SIMBAD_$precision" };
    }

    push @cat, $self->get_sdss();

    # Add Y,Z angle for database
    my $cat_id = 1;
    foreach $cat (@cat) {
	$cat->{id} = $cat_id++;
	$cat->{obsid} = $self->obsid;
	($cat->{y_angle}, $cat->{z_angle}) = main::calc_y_z_r($self->ra, $self->dec, $self->roll,
							      $cat->{ra}, $cat->{dec});
    }

    return \@cat;
}

##************************************************************************---
sub get_sdss {
#
#
# SELECT top 20 p.objid, p.run, p.rerun, p.camcol, p.field, p.obj, p.type, 
#        p.ra, p.dec, p.u,p.g,p.r,p.i,
#        p.z, p.Err_u, p.Err_g, p.Err_r,p.Err_i,p.Err_z 
#  FROM fGetNearbyObjEq(180,-0.2,3) n, PhotoPrimary p 
#  WHERE n.objID=p.objID
##************************************************************************---
    my $self = shift;
    my @cat;
    my %sdss_opt = (format => 'xml',
		    topnum => 300,
		    ra     => sprintf("%.5f",$self->ra),
		    dec    => sprintf("%.4f",$self->dec),
		    radius => $self->extr_rad,
		   );
    my $url = $par{sdss_url}
      . '/x_radial.asp?'
      . join('&', map { "$_=$sdss_opt{$_}" } keys %sdss_opt);
    $log->message("Getting SDSS stellar objects: $url");
    my ($xml, $error) = get_url($url, timeout => 120);
    die "Error accessing '$url': $error\n" if defined $error;

    my $vals = XMLin($xml);
    return () unless ref($vals->{Answer}{Row}) eq 'ARRAY';

    foreach my $object (@{$vals->{Answer}{Row}}) {
	next unless $object->{type} == 6; # Accept only stellar objects
	my $ra = $object->{ra};
	my $dec = $object->{dec};
	my ($rah, $decd) = dec2hms($ra, $dec);
	push @cat, {name    => $object->{objid},
		    ra_hms  => $rah,
		    dec_hms => $decd,
		    ra      => $ra,
		    dec     => $dec,
		    mag     => $object->{r} ,
		    catalog => "SDSS" };
    }

    return @cat;
}


##************************************************************************---
sub clean_up {
##***************************************************************************
    my $self = shift;
    (my $clean = shift) =~ tr/A-Z/a-z/;

    my $asp_glob = "ASP_L1_STD* aca0* asp05* pcad0* sim05* obspar*"; # pcad*asol*";
    my %glob_list = (event => "*evt2.fits* $asp_glob",
		     all   => "*src2.fits* *evt2.fits* axaf*obs0a.par* $asp_glob",
		     most  => "hrc*evt2.fits* acis*evt2.fits* $asp_glob",
		    );
    my $glob = $glob_list{$clean} || '';
	
    foreach (glob($glob)) {
	$log->message("Cleaning up $_");
	my $io = io($_);
	$io->type eq 'dir' ? $io->rmtree : $io->unlink;
    }
}

##************************************************************************---
sub gzip_fits {
##***************************************************************************
    my $self = shift;

    for (glob("*.fits")) {
	system("gzip $_") and die "Couldn't gzip $_" ;
    }
}

##*************************************************************************--
sub preclean {
# Preclean files.  Useful to force cross correlation or new catalog
# for already processed obsid
##***************************************************************************
    my $self = shift;

    if (my $file_list = shift) {
	foreach my $file (glob($file_list)) {
	    if (-e $file) {	# remember glob('doesntexist') returns 'doesntexist'
		$log->message("Removing $_");
		unlink $_ ;
	    }
	}
    }
}

##************************************************************************---
sub search_hash {
##***************************************************************************
#    my ($table, $field_values) = @_;
#    my @fields = sort keys %$field_values; # sort required
#    my @values = @{$field_values}{@fields};
#    my $qualifier = "";
#    $qualifier = "where ".join(" and ", map { "$_=?" } @fields) if @fields;
#    $sth = $dbh->prepare_cached("SELECT * FROM $table $qualifier");
#    return $dbh->selectall_arrayref($sth, {}, @values);
}


##****************************************************************************
##****************************************************************************
package XraySource;
##****************************************************************************
##****************************************************************************

use Quat;
use Ska::Convert qw(dec2hms);
use Data::Dumper;

# A pseudo-object to provide a clean way of getting all the fields
# needed for the astromon_xray_source table

use Class::MakeMethods::Standard::Hash ( 
					scalar => [ qw(ra dec net_counts snr double_id) ],
					array => 'xray_sources',
					hash => 'src2',
					object => 'obs',
				       );

# 			   obsid      => 'int',
# 			   id         => 'varchar(22)', # ala 'CXOU J123456.7+765432'
# 			   ra         => 'double precision', 
# 			   dec        => 'double precision',
# 			   net_counts => 'float',  
# 			   net_rate   => 'float',  
# 			   y_angle   => 'float',  
# 			   z_angle   => 'float',  
# 			   r_angle   => 'float',  
# 			   snr        => 'float',
# 			   status_id  => 'int null',
#

##************************************************************************---
sub new {
##************************************************************************---
    my $class = shift;
    my $self = { @_ };
    bless $self, $class;

    $self->$_( $self->src2($_) ) for qw(ra dec net_counts snr);
    
    return $self;
}


##************************************************************************---
sub obsid {
##***************************************************************************
    my $self = shift;
    return $self->{obsid} = $self->obs->obsid;
}

##************************************************************************---
sub name {
##***************************************************************************
    my $self = shift;
    local $_;
    my ($ra_hms, $dec_hms) = dec2hms($self->src2('ra'), $self->src2('dec'));
    $ra_hms =~ s/(\..).*/$1/;
    $dec_hms =~ s/\..*//;
    (my $name = "CXOU ${ra_hms}${dec_hms}") =~ s/://g;
    return $self->{name} = $name;
}
    
##************************************************************************---
sub y_angle {
##***************************************************************************
    my $self = shift;
    return $self->{y_angle} if defined $self->{y_angle};

    my $obs = $self->obs;
    ($self->{y_angle}, $self->{z_angle}, $self->{off_axis_angle}) =
      main::calc_y_z_r($obs->ra, $obs->dec, $obs->roll, $self->src2('ra'), $self->src2('dec'));
    
    return $self->{y_angle};
}

##************************************************************************---
sub z_angle {
##***************************************************************************
    my $self = shift;
    return $self->{z_angle} if defined $self->{z_angle};

    my $obs = $self->obs;
    ($self->{z_angle}, $self->{z_angle}, $self->{r_angle}) =
      main::calc_y_z_r($obs->ra, $obs->dec, $obs->roll, $self->src2('ra'), $self->src2('dec'));
    
    return $self->{z_angle};
}

##************************************************************************---
sub r_angle {
##***************************************************************************
    my $self = shift;
    return $self->{r_angle} if defined $self->{r_angle};

    my $obs = $self->obs;
    ($self->{r_angle}, $self->{z_angle}, $self->{r_angle}) =
      main::calc_y_z_r($obs->ra, $obs->dec, $obs->roll, $self->src2('ra'), $self->src2('dec'));
    
    return $self->{r_angle};
}

##************************************************************************---
sub id {
##***************************************************************************
    my $self = shift;
    return $self->src2('component');
}

##************************************************************************---
sub status_id {
##***************************************************************************
    my $self = shift;
    return $self->{status_id};
}

##************************************************************************---
sub near_neighbor_dist {
##***************************************************************************
    my $self = shift;
    
    return $self->{near_neighbor_dist} if defined $self->{near_neighbor_dist};

    # Find the distance to the nearest neighbor
    my $near_neighbor_dist2 = 1e38;
    foreach my $xs (@{$self->xray_sources}) {
	next if $xs == $self;
	if ((my $dist2 = ($xs->y_angle - $self->y_angle)**2 + ($xs->z_angle - $self->z_angle)**2) 
	    < $near_neighbor_dist2) {
	    $near_neighbor_dist2 = $dist2;
	}
    }

    return $self->{near_neighbor_dist} =  sqrt($near_neighbor_dist2) ;
}



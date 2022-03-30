package CalcOffset;

##***********************************************************************
#
# calc_offset.pl - Calculate offsets, based on fix_offset.pl
#
# Copyright 2001  Tom Aldcroft
#
##***********************************************************************

use warnings;
use Text::ParseWords;
use Expect::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::Slatec;
use App::Env;
use Cwd;
use CFITSIO::Simple;
use Data::Dumper;

1;

##****************************************************************************
# MAIN PROGRAM
##****************************************************************************

sub calc_offset {
    %par = (loud => 1,
	    archive => 1,
	    @_);

    $work_dir = cwd;

    check_options();
    get_obspar();
    get_current_calalign();
    get_acal_fidpr();
    return calculate_offset();
}

##****************************************************************************
sub check_options {
##****************************************************************************
    local $_;

    $par{obsid} = $par{obs_id} + 0;
    ($yr, $mon, $day) = ($1,$2,$3) if ($par{'date-obs'} =~ /(\d+)-(\d+)-(\d+)/);
    $val{RA} = $par{ra_nom};
    $val{DEC} = $par{dec_nom};

    fatal_error("Header information does not contain DATE keyword")
	unless (defined $yr && defined $mon && defined $day);
    fatal_error("Processing date is $yr-$mon-$day, but observations processed \n"
		."before 2000-Feb-01 cannot be fixed by this tool.  Reprocessing is\n"
		."required to remove any offset. Please contact the CXC help desk\n"
		."or your contact scientist to arrange for reprocessing")
	if ($yr < 2000 || ($yr == 2000 and $mon < 2));

    fatal_error("Header information does not contain OBS_ID keyword")
	unless ($par{obsid});
    fatal_error("Header information does not contain ACSDSVER keyword")
	unless ($par{ascdsver});
    $par{ascdsver} =~ s/\s*$//;	# Get rid of any trailing spaces
    fatal_error("Invalid obsid: $par{obsid}") unless ($par{obsid} =~ /\d+/);

    fatal_error("Couldn't parse RA and DEC from header information\n\n$par{header}") 
	unless (defined $val{RA} and defined $val{DEC});

    # Remove leading or trailing spaces
    foreach (keys %par) {
	$par{$_} =~ s/\A\s+|\s+\Z//g if defined $par{$_};
    }
    print "\n" if ($par{loud});

    return 1;
}

##****************************************************************************
sub get_obspar {
##****************************************************************************
# Get obspar file, either from current directory or from archive

    (@obspars) = get_archive_files(obsid     => $par{obsid},
				   prod      => "obspar",
				   file_glob => "axaf*[0-9]_obs0a.par*",
				   version   => ['all'],
				   );
    fatal_error("No obspar files found") unless (@obspars);

# Parse obspar

    foreach (@obspars) {
	%obspar = read_param_file ($_);
	last if ($obspar{ascdsver} eq $par{ascdsver});
    }

    fatal_error("Obsid mismatch $par{obsid} != $obspar{obs_id}") unless ($par{obsid} == $obspar{obs_id});

#     ($detector = $obspar{detector}) =~ tr/A-Z/a-z/;

    return 1;
}

##****************************************************************************
sub get_acal_fidpr {
##****************************************************************************
    my $using_most_recent = 0;

    (@fidpr)= get_archive_files(obsid     => $par{obsid},
				prod      => "asp1[*fidpr*]",
				file_glob => "pcad*fidpr*fits*",
				version   => ['all'],
				ascdsver  => $par{ascdsver},
				);

    # Check if we got any fidprops files with the specified version
    unless (@fidpr) {

	# No luck, so try to get ANY files and choose the highest version

	(@fidpr)= get_archive_files(obsid     => $par{obsid},
				    prod      => "asp1[*fidpr*]",
				    file_glob => "pcad*fidpr*fits*",
				    version   => ['all'],
				    );

	# Still can't find anything?
	fatal_error("No FIDPROPS files found for this Obsid") unless (@fidpr);
	    
	$using_most_recent = 1;

	# Pick out the files with the highest version number
	my $version = 0;
	map { $version = $1 if (/pcad\S*\d+N(\d+)/ && $1 > $version) } @fidpr;
	@fidpr = grep { /pcad\S*\d+N(\d+)/ && $1 == $version } @fidpr;
    }

    # Read in the header keywords
    ($fid_hdr, $error) = fits_read_hdr($fidpr[0], 'FIDPROPS');
    fatal_error("Could not read FIDPROPS header from $fidpr[0]") if ($error);
    $sim_z_offset = $fid_hdr->{LSI0STT3} + $fid_hdr->{STT0STF3};

    # If processing fell through to using the most recent version, then print
    # a big warning message and set $par{ascdsver} to the FIDPR ASCDSVER

    if ($using_most_recent) {
	$msg = <<MSG_END
No FIDPROPS files with version='$par{ascdsver}'.  
This may be because you processed the event file with a CIAO tool which changed 
the version header keyword (ASCDSVER), or it may be because the event file was
reprocessed without the aspect being reprocessed.

IMPORTANT:
 Fix_offset will now use the most recent available aspect files:
  Processing date    : $fid_hdr->{DATE}
  Processing version : $fid_hdr->{ASCDSVER}
  
This assumes that the event file is the latest available version
currently in the archive, and that none of the coordinate reference
keywords in the file have been changed.  If this is not the case, the
resultant dmhedit commands may be incorrect.

MSG_END
    ;
	print $msg;
	($par{ascdsver} = $fid_hdr->{ASCDSVER}) =~ s/\A\s+|\s+\Z//g;;
    }

    # Make sure RA_NOM and DEC_NOM are defined

    fatal_error("FIDPROPS file did not contain nominal RA and Dec")
	unless (defined $fid_hdr->{RA_NOM} and defined $fid_hdr->{DEC_NOM});

    # Make sure that RA_NOM and DEC_NOM match RA and DEC in the event file

    if (3600*radec_dist($val{RA}, $val{DEC}, $fid_hdr->{RA_NOM}, $fid_hdr->{DEC_NOM}) > .01) {
	$msg = <<ERROR_END
The RA_NOM and DEC_NOM in the archive aspect fid properties do not match the
values found in the event file header.  This means the event file coordinate
reference has been changed since retrieval from the archive.  This could happen
if you ran reproject_events, or else from running this thread.  It is important
for this tool that the coordinate reference keywords be unchanged from those
created in pipeline processing.  

Note that if you have already run this thread and fixed your event file, 
you cannot re-run it to see if it was fixed properly.  

The current values are: 
    Event RA, DEC = $val{RA}, $val{DEC}
    Fidprops RA, DEC = $fid_hdr->{RA_NOM}, $fid_hdr->{DEC_NOM}
ERROR_END
    ;
	fatal_error($msg);
    }

    # Find the last valid fidprops file based on a sort index defined
    # as <version><time>, so that a later version is always a later time.
    # Define the sort key, Sort the (0 .. $#fid_key) array by the sort key,
    # and find the last element

    my @fid_key = map { /pcad\S*(\d+)N(\d+)/; "$2$1" } @fidpr;
    @fid_index = sort { $fid_key[$a] <=> $fid_key[$b] } (0 .. $#fid_key);
    $fidpr = @fidpr[$fid_index[-1]];
    
    # Read the pos_lsi, fid_num, and fid_status columns from that file
    ($fid_pos_lsi1, $fid_num1, $fid_status1) = fits_read_bintbl($fidpr, qw(P_LSI ID_NUM ID_STATUS),
								{ extname => 'FIDPROPS', rethash => 0 });
    fatal_error("Could not read data from fid properties file $fidpr")
      unless (defined $fid_pos_lsi1);
    
    (@acal)= get_archive_files(obsid     => $par{obsid},
			       prod      => "asp1[*acal*]",
			       file_glob => "pcad*acal*fits*",
			       version   => ['all'],
			       ascdsver  => $par{ascdsver},
			       );

    unless (@acal) {
	$msgs = <<ERROR_END
No ACACAL files with version='$par{ascdsver}'.
ERROR_END
    ;
	fatal_error($msgs);
    }

    # for ACA CAL, it is sufficient to read the first file

    ($aca_misalign) =
	fits_read_bintbl($acal[0],('ACA_MISALIGN'), { extname => 'ACACAL', rethash => 0 });
    fatal_error("Could not read data from ACA CAL file $acal[0]")
	unless (defined $aca_misalign);
    $aca_misalign1 = $aca_misalign(;-);
}

##****************************************************************************
sub get_current_calalign {
##****************************************************************************
    return () unless ($obspar{"date-obs"} =~ /(.+)T(.+)/);
    my ($date, $time) = ($1,$2);

    # Get the latest alignment for this time and date, in format "file <number>"

    my $ascds_env = App::Env->new('ASCDS');
    local %ENV = %{$ascds_env};
    $ENV{PFILES} = "$ENV{PWD};$ENV{PFILES}";

    if ($par{caldb}) {
	$ENV{CALDB} = $par{caldb};
	$ENV{CALDBCONFIG} = "$par{caldb}/software/tools/caldb.config";
	$ENV{CALDBALIAS} = "$par{caldb}/software/tools/alias_config.fits";
    }

    $_ = `quzcif CHANDRA PCAD - - align $date $time -`;

    my ($align) = split;
    fatal_error("No CALALIGN file found for date $date and time $time") unless (-e $align);

    my ($instr_id, $aca_misalign) =
	fits_read_bintbl($align, ('INSTR_ID','ACA_MISALIGN'), { extname => 'CALALIGN',
								rethash => 0 });

    fatal_error("Could not read alignment data from CALALIGN file $align")
	if ($status);

    $i_det = 0;
    while ($instr_id->[$i_det] ne $obspar{detector} and $i_det < 4) {
	$i_det++;
    }

#    print "align = '$align'\n";
#    print "aca_misalign = $aca_misalign\n";
#    print "idet = $i_det, dims aca_misalign = ", $aca_misalign->dims(), "\n";

    $aca_misalign2 = $aca_misalign(:,:,$i_det;-);

    # Read fid light positions from second extension of CALALIGN file

    ($fid_num2, $fid_si2, $fid_pos_lsi2, $fid_corr2{y}, $fid_corr2{z}, $fid_lim2{y}, $fid_lim2{z}) =
	fits_read_bintbl($align, ('FID_NUM','FID_SI','FID_POS_LSI', 'FID_Y_CORR', 'FID_Z_CORR',
				  'FID_Y_LIM', 'FID_Z_LIM'),
				  { extname => 'CALALIGN1',
				    rethash => 0 });

    fatal_error("Could not read fid data from CALALIGN file $align")
	if ($status);

}

##****************************************************************************
sub get_archive_files {
##****************************************************************************
    my %arg = (version => ['last'],
	       dir     => $work_dir,
	       loud    => $par{loud},
	       @_);
    local $_;

    my @files;

    # Check if file is already available

    unless ((@files = glob($arg{file_glob})) or !$par{archive}) {

	# Start arc5gl if necessary
	unless ($arc5gl) {
	    print "Starting arc5gl\n" if ($par{loud});
	    $arc5gl = new Expect::Simple { Cmd => "arc5gl -stdin",
					   Prompt => 'ARC5GL> ',
					   DisconnectCmd => 'exit',
					   Verbose => 0,
					   Debug => 0,
					   Timeout => 240
					   };
	}

	print "Getting $arg{prod} from archive\n";
	$arc5gl->send("loud");
	$arc5gl->send("cd $arg{dir}") if ($arg{dir});

	if ($arg{obsid}) {
	    $arc5gl->send("obsid=$arg{obsid}");
	} elsif ($arg{tstart}) {
	    $arc5gl->send("tstart=$arg{tstart}");
	    $arc5gl->send("tstop=$arg{tstop}");
	} else {
	    fatal_error("get_archive_files:: Need obsid or tstart");
	}
	
	for $version (@{$arg{version}}) {
	    $arc5gl->send("version=$version");
	    $arc5gl->send("get $arg{prod}");
	    @files = glob($arg{file_glob});
	    last if (@files);
	}
    }
	
    # Check the ASCDS version if needed

    if (defined $arg{ascdsver}) {
	my @files_good = ();
	foreach (@files) {
	    my ($keys, $status) = fits_read_hdr($_,2);
	    next if ($status or not defined $keys->{ASCDSVER});
	    $keys->{ASCDSVER} =~ s/\A\s+|\s+\Z//g;
	    push @files_good, $_ if ($keys->{ASCDSVER} eq $arg{ascdsver});
	}
	@files = @files_good;
    }

    return @files;
}

##****************************************************************************
sub calculate_offset {
##****************************************************************************
    $r2a = 206264.8;
    $d2r = 3.14159265 / 180.;
    $align_diff = transpose($aca_misalign2) x $aca_misalign1 * $r2a;

    for $i (0 .. $#{$fid_status1}) {
	next unless ($fid_status1->[$i] eq 'GOOD');
	push @y1, $fid_pos_lsi1->at(1,$i);
	push @z1, $fid_pos_lsi1->at(2,$i);
	$n_fid1++;
	for $j (0 .. nelem($fid_num2)-1) {
	    next unless ($fid_si2->[$j] eq $obspar{detector}
			 && $fid_num2->at($j) == $fid_num1->at($i));
	    $n_fid2++;

	    # Calculate the fid position correction term
	    for $dir qw(y z) {
		$fid_corr{$dir} = 0;
		map { $fid_corr{$dir} +=  $fid_corr2{$dir}->at($_,$j) * $sim_z_offset**$_ } (0..4);
		$fid_corr{$dir} = $fid_lim2{$dir}->at(0,$j) if ($fid_corr{$dir} < $fid_lim2{$dir}->at(0,$j));
		$fid_corr{$dir} = $fid_lim2{$dir}->at(1,$j) if ($fid_corr{$dir} > $fid_lim2{$dir}->at(1,$j));
	    }

	    #	    print "fid corr{y} = $fid_corr{y} fid corr{z} = $fid_corr{z} \n";
	    #	    print("fid pos_lsi2 for fid ",$fid_num2->at($j)," = ",$fid_pos_lsi2->at(1,$j)-$fid_corr{y},", ", 
	    #	          $fid_pos_lsi2->at(2,$j)-$fid_corr{z}, "\n");

	    push @y2, $fid_pos_lsi2->at(1,$j) - $fid_corr{'y'};
	    push @z2, $fid_pos_lsi2->at(2,$j) - $fid_corr{'z'};
	}
    }

    fatal_error("Could not find all fid lights in CALALIGN file")
      unless ($n_fid1 == $n_fid2 && $n_fid1 != 0);

    # Try to calculate the differences using the "correct" approach
    # (analytical least squares fit method)

    ($y_diff, $z_diff) = calc_fid_dy_dz(\@y1, \@z1, \@y2, \@z2);

    # If it didn't work, just take average of diffs then convert to arcsec

    unless (defined $y_diff and defined $z_diff) {
	$y_diff = (sum(pdl @y1) - sum(pdl @y2)) / $n_fid1 * 20;  
	$z_diff = (sum(pdl @z1) - sum(pdl @z2)) / $n_fid1 * 20;  
    }

    $dy = $align_diff->at(1,0) - $y_diff;
    $dz = $align_diff->at(2,0) - $z_diff;
	    
    $rollr = $obspar{roll_nom} * $d2r;
    $decr  = $val{DEC} * $d2r;

    if (abs($val{DEC}) < 89) {
	$d_RA = ($dy * cos($rollr) - $dz * sin($rollr)) / cos($decr) / 3600.;
	$d_dec= ($dy * sin($rollr) + $dz * cos($rollr)) / 3600.;
    }

    return ($dy, $dz, $d_RA, $d_dec);
}


##***************************************************************************
sub calc_fid_dy_dz {
##  Calculate dy and dz offset due to change in fid coordinates
##***************************************************************************
    my ($y1, $z1, $y2, $z2) = @_;

    my $n_fid = @{$y1};
    my $Mi = zeroes(3,3);
    my $Xi = zeroes(3);

#   print "@{$y1}, @{$z1}, @{$y2}, @{$z2}\n";

    return () unless ($n_fid >= 2);

    for (0 .. $n_fid-1) {
	my $y_obs = $y1->[$_];
	my $z_obs = $z1->[$_];

	$y = $y2->[$_];
	$z = $z2->[$_];

	$Mi = $Mi + pdl [[ 1, 0.0, -$z ], 
			 [ 0.0, 1, $y  ], 
			 [-$z, $y, ($y**2+$z**2)]];

	$Xi = $Xi -  pdl [ $y - $y_obs ,  $z - $z_obs, $z * $y_obs -   $y * $z_obs];
    }

    $delta =  matinv($Mi) x transpose($Xi);
#    print "Mi = $Mi\nXi = $Xi\ndelta = $delta\n";
    
    return ($delta->at(0,0)*20, $delta->at(0,1)*20);  # Return dy and dz in arcsec
}

##****************************************************************************
sub fatal_error {
##****************************************************************************
    my $msg = shift;

    die "ERROR: $msg\n";
}

###################################################################################
sub read_param_file {
###################################################################################
    my $file = shift;
    my %param = ();
    $file = "gunzip --stdout $file |" if ($file =~ /(\.gz|\.Z|\.z)$/);
    open (PAR, $file) || die "Couldn't open parameter file '$file'\n";
    while (<PAR>) {
	@fields = quotewords(",", 0, $_);
	$param{$fields[0]} = $fields[3];
    }
    close PAR;

    return %param;
}

##**********************************************************************
sub radec_dist {
##**********************************************************************
    my ($a1, $d1, $a2, $d2) = @_;
    my $d2r = 3.14159265358979/180.;

    return(0.0) if ($a1==$a2 && $d1==$d2);

    return acos( cos($d1*$d2r)*cos($d2*$d2r) * cos(($a1-$a2)*$d2r) +
              sin($d1*$d2r)*sin($d2*$d2r)) / $d2r;
}

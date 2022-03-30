use Ska::RunAsp;

my %ascds;

##****************************************************************************
sub reprocess_aspect {
##****************************************************************************
    my $asol;

    unless (%ascds) {
	%ascds = grabenv('tcsh', 'source /home/ascds/.ascrc -r release');
	$ascds{ENVIRONMENT_NAME} = 'ascds';
    }

    if ($par{reprocess}) {
	my $range;
	local %ENV = %ascds;
	print STDERR "Reprocessing aspect -- this may take awhile...\n" if ($par{loud});

	# pre-clean if needed, and set parameter file to defaults
	system("rm -rf ASP_L1_STD") if (-e "ASP_L1_STD");
	print STDERR "Using default ASCDS asp_l1_std.par\n";
	system("punlearn asp_l1_std");
	$ENV{PWD} = $work_dir;

	# allow for a custom CALDB in the current work directory
	if (-e "$HOME/CALDB") {  
	    $ENV{CALDB} = "$HOME/CALDB";
	    print "Using local CALDB: $HOME/CALDB\n";
	}

	# Make sure that range doesn't exceed obs length
	$range = ($par{range} and $par{range} < $obspar{tstop} - $obspar{tstart}) ?
	    "0:+$par{range}" : "0:";

	# Run aspect pipeline
	&Ska::RunAsp::go(obsid  => $obsid,
		    clean  => 0,
		    range  => $range,
		    qual   => "obs${obsid}_qual.ps",
		    param  => $par{repro_param},
		    );

	# Get rid of existing fidpr and acal, and copy aspect solution, 
	# fidpr, and acal into work dir

	chdir $work_dir;
	$ENV{PWD} = $work_dir;
	unlink glob "pcad*asol1.fits* pcad*fidpr1.fits* pcad*acal1.fits*";

	system("cp ASP_L1_STD/out1/pcad*asol1.fits ./");
	system("cp ASP_L1_STD/out1/pcad*fidpr1.fits ./");
	system("cp ASP_L1_STD/out1/pcad*acal1.fits ./");
    }

    # Find aspect solution for reproject_events
    if ($par{reprocess} || $par{reapply}) {
	($asol) = glob("pcad*asol1.fits ASP_L1_STD/out1/pcad*asol1.fits");
	return () unless ($asol);
	print "Using aspect solution file $asol\n" if ($par{loud});
    }

    return $asol;
}


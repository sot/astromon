# Original calling sequence in celmon
#    # Optionally make gif images of each source match using ds9
#    make_images(@matches) if ($par{make_images} and $evt2 and -e $evt2); 

##****************************************************************************
sub make_images {
##****************************************************************************
    eval 'use Image::DS99';
    my @sources = @_;
    local $_;
    local %ENV = %ENV_plain;

    my $good_pos_ref = qr/tycho|simbad_high|icrs/i;
    return 1 unless (grep {$_->{pos_ref} =~ $good_pos_ref} @sources);

    # Define various file names for use later
    my $ds9_pid;

    return 1 if $ds9_failed;	# Already tried and failed to start ds9

    my $ds9 = Image::DS9->new( { Server => 'Celmon_DS9' });
    unless ( $ds9->nservers ) {
	message("Starting ds9...\n");
	my $geometry = "617x728-0+0";
	$ds9_pid = open DS9, "ds9 -geometry $geometry -title Celmon_DS9 2>&1 |";
	$ds9->wait(30) or do {
	    message("ERROR - Could not start ds9\n");
	    $ds9_failed = 1;
	    return 1;		# Can still do other processing without ds9
	}
    }

    $ds9->zoom(to => 1); 
    my $file = (-e "source_evt2.fits") ? "source_evt2.fits" : $evt2;
    $ds9->file($file);

    my $n_source = 0;
    foreach $source (@sources) {
	$n_source++;
	next unless $source->{pos_ref} =~ $good_pos_ref;
	my ($ra, $dec) = ($source->{cat_ra}, $source->{cat_dec});
	my $psfile = "$work_dir/img_${n_source}.ps";
	my $sclfac = ($source->{detector} =~ /hrc/i) ? 4 : 1;

	my $tmp_reg = tmpnam();
	my @reg;
	push @reg, "# Region file format: DS9 version 3.0\n";
	push @reg, "global color=green\n";
	push @reg, sprintf('fk5;circle(%.6f,%.6f,1.0")%s', $ra, $dec,"\n");
	push @reg, sprintf('fk5;box(%.6f,%.6f,0.5",0.5",45)%s', $source->{ra}, $source->{dec},"\n");
	write_file($tmp_reg, @reg);

	$ds9->zoom(to => 1); 
	$ds9->bin(factor => 16);
	$ds9->pan(to => ($ra,$dec), "wcs", "fk5");
	$ds9->bin(factor => 0.5 * $sclfac);
	$ds9->zoom(to => 8);
	$ds9->scale("log");
	$ds9->scale(mode => "minmax");
	$ds9->cmap("Heat");
	$ds9->print(destination => 'file');
	$ds9->print(filename => $psfile);
	$ds9->print(resolution => 75);
	$ds9->regions("deleteall");
	$ds9->regions(load => $tmp_reg);
	$ds9->print();
	unlink $tmp_reg;
	system("ps2gif $psfile");
	unlink $psfile;
    }

    # Close down ds9 each time
    kill 9 => $ds9_pid;
    message("Killed ds9 pid=$ds9_pid\n");
    sleep 4;
    close DS9;

    1;
}


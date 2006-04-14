##****************************************************************************
sub filter_source_list {
##****************************************************************************
    my $self = shift;
    
    # Read the cell detect source file, and then eliminate any sources
    # which are within $par{excl_rad} of another source

    my @src_good = ();
    my $src_data = parse_table($self->src2); # Read src2 FITS file as array of hashes
    foreach my $src1 (@{$src_data}) {
	my $bad = 0;
	foreach my $src2 (@{$src_data}) {
	    next if $src1 == $src2;
	    $bad = 1 if (radec_dist($src1->{ra}, $src1->{dec}, $src2->{ra}, $src2->{dec})*3600 < $par{excl_rad});
	}
	push @src_good, $src1 unless $bad;
    }
    $self->sources( \@src_good );
}


##************************************************************************---
sub get_cell_detect {
##***************************************************************************
    my $self = shift;
# Read the cell detect source file, and then eliminate any sources
# which are within $par{excl_rad} of another source
    my $excl_rad = shift;

    my @celldet_good = ();
    $cmd = q{dmlist "} . $self->src2 . q{[col ra,dec,net_counts,snr]" data,clean};
    $log->message("$cmd");
    my @celldet = `$cmd`;
    @celldet = grep !/^\#/, @celldet;
    for my $i (0 .. $#celldet) {
	my $bad = 0;
	for my $j (0 .. $#celldet) {
	    next if ($i == $j);
	    my ($ra1, $dec1) = split ' ', $celldet[$i];
	    my ($ra2, $dec2) = split ' ', $celldet[$j];
	    $bad = 1 if (radec_dist($ra1, $dec1, $ra2, $dec2)*3600 < $excl_rad);
	}
	push @celldet_good, $celldet[$i] unless ($bad);
    }
    
    return @celldet_good;
}

##***************************************************************************
sub xcorr_sources {
##***************************************************************************
    my $self = shift;

    my $n_match = 0;
    my @matches = ();
    my ($offset_dy, $offset_dz, $offset_RA, $offset_dec) = (0,0,0,0);
    
    my @celldet = get_cell_detect($par{excl_rad});

    # Calculate offset for specified aspect products based on date and obsid
    #  THIS CODE DISABLED FOR NOW
    if (0) {
	($offset_dy, $offset_dz, $offset_RA, $offset_dec) = CalcOffset::calc_offset(%obspar,
										    caldb => $par{caldb});
	printf "Offsets dy,dz = %.2f,%.2f RA,dec= %.5f,%.5f\n", $offset_dy, $offset_dz,
	  $offset_RA, $offset_dec;
    }

  SOURCE: foreach (@celldet) {
	my ($ra, $dec, $counts, $snr) = split;
	$ra += $offset_RA;
	$dec += $offset_dec;
      
      CAT: foreach my $cat (@cat) {
	    my $delta_rad = 3600 * radec_dist($ra, $dec, $cat->{ra}, $cat->{dec});
	    next CAT unless ($delta_rad < $par{det_rad});
	
	    my $delta_dec = ($cat->{dec} - $dec) * $d2a;
	    my $delta_ra  = ($cat->{ra} - $ra) * cos($cat->{dec}*$d2r) * $d2a;

	    my $cr = cos($obspar{roll_nom} * $d2r);
	    my $sr = sin($obspar{roll_nom} * $d2r);
	    my $dz = $cr * $delta_dec - $sr * $delta_ra; #  - $offset_dz;   
	    my $dy = $sr * $delta_dec + $cr * $delta_ra; #  - $offset_dy;

	    my $delta_field_dec = ($obspar{dec_nom} - $dec) * $d2m;
	    my $delta_field_ra  = ($obspar{ra_nom} - $ra) * cos($obspar{dec_nom}*$d2r) * $d2m;
	    my $field_dz = $cr * $delta_field_dec - $sr * $delta_field_ra;   
	    my $field_dy = $sr * $delta_field_dec + $cr * $delta_field_ra;
	    my $field_dr = sqrt($field_dy**2 + $field_dz**2);

	    my $cnt_rate = $counts / (($obspar{tstop} - $obspar{tstart}) / $self->frame_time);

	    if ($obspar{ascdsver} =~ /R4CU(\d)UPD([\d\.]+)/) {
		$version = sprintf "%.3f", $1 + $2/100.0 ;
	    } elsif ($obspar{ascdsver} =~ /(\d+)\.(\d+)\.(\d+)/) {
		$version = sprintf "%.4f", $1 + $2/100.0 + $3/10000.0;
	    } else {
		$version = "0.000";
	    }

	    my %out = ();
	    $out{obsid}    = sprintf("%d", $obsid);
	    $out{detector} = $obspar{detector};
	    $out{target}   = $obspar{object};
	    $out{object}   = $cat->{name};
	    $out{date_obs} = $obspar{'date-obs'};
	    $out{grating}  = $obspar{'grating'};
	    $out{sim_z}    = sprintf("%.2f", $obspar{sim_z});
	    $out{ra}       = sprintf("%.7f", $ra);
	    $out{dec}      = sprintf("%.7f", $dec);
	    $out{cat_ra}   = sprintf("%.7f", $cat->{ra});
	    $out{cat_dec}  = sprintf("%.7f", $cat->{dec});
	    $out{dy}       = sprintf("%.2f", $dy);
	    $out{dz}       = sprintf("%.2f", $dz);
	    $out{dr}       = sprintf("%.2f", $delta_rad);
	    $out{counts}   = sprintf("%d",   $counts);
	    $out{cnt_rate} = sprintf("%.1e", $cnt_rate);
	    $out{field_dy} = sprintf("%.2f", $field_dy);
	    $out{field_dz} = sprintf("%.2f", $field_dz);
	    $out{field_dr} = sprintf("%.2f", $field_dr);
	    $out{snr}      = sprintf("%.1f", $snr);
	    $out{tstart}   = sprintf("%d",   $obspar{tstart});
	    $out{fids}     = $self->fids;
	    $out{ascdsver} = $obspar{ascdsver};
	    $out{version}  = $version;
	    $out{pos_ref}  = $cat->{catalog};
	    $out{status}   = '';

	    push @matches, { %out };

	    $n_match++;	
	    $match = sprintf("%5d %-10s %-10s %-14s %12.7f %12.7f %6.2f %6.2f %7.1e %6.1f %6.2f %6.2f",
 			     $obspar{obs_id}, $detector, $cat->{catalog}, $cat->{name}, $ra, $dec,
			     $dy, $dz, $cnt_rate, $snr, $field_dy, $field_dz);
	    $log->message($match);
	} 
    }  

    $log->message("OBSID " . $self->obsid . ": Successful cross-correlation for $obspar{object} ($n_match matches)");
    $self->n_match($n_match);

    return @matches;
}


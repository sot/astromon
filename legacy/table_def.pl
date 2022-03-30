# -*- cperl -*-
    (
     astromon_obs => [		# observation
		      obsid      => 'int',
		      version    => 'float null',
		      detector   => 'varchar(6) null',  
		      target     => 'varchar(28) null', 
		      grating    => 'varchar(4) null',  
		      sim_z      => 'float null',  
		      date_obs   => 'varchar(20) null', # 2006-04-13T20:23:57
		      tstart     => 'double precision null',  
		      fids       => 'varchar(20) null', 
		      ascdsver   => 'varchar(32) null', 
		      ra_nom     => 'double precision null',
		      dec_nom    => 'double precision null',
		      roll_nom   => 'double precision null',
		      process_status => 'varchar(128) null',
		     ],

     # detected x-ray source 
     astromon_xray_src => [
			   obsid      => 'int',
			   version    => 'float null',
			   id         => 'varchar(22)', # ala 'CXOU J123456.7+765432'
			   ra         => 'double precision', 
			   dec        => 'double precision',
			   skyx       => 'float',       
			   skyy       => 'float',       
			   counts     => 'float',  
			   cnt_rate   => 'float',  
			   field_dy   => 'float',  
			   field_dz   => 'float',  
			   field_dr   => 'float',  
			   snr        => 'float',
			   status_id  => 'int null',
			  ],
     # catalog source
     astromon_cat_src => [
			  obsid      => 'int',
			  catalog    => 'varchar(16)', 
			  id         => 'varchar(24)',
			  ra         => 'double precision', 
			  dec        => 'double precision',
			  mag        => 'float null',
			  skyx       => 'float',
			  skyy	  => 'float',
			 ],

     # Cross-correlated source
     astromon_xcorr => [
			obsid   => 'int',
			catalog => 'varchar(16)',
			c_id    => 'varchar(24)',
			x_id    => 'int',
			dy      => 'float',  
			dz      => 'float',  
			dr      => 'float',
		       ],
     # Status messages
     astromon_status => [
			 id     => 'int',
			 status => 'varchar(256)',
			],
    );


CREATE TABLE astromon_obs (
    obsid      int not null,
    version    float null,
    detector   varchar(6) null,  
    target     varchar(28) null, 
    grating    varchar(4) null,  
    sim_z      float null,  
    date_obs   varchar(20) null, -- 2006-04-13T20:23:57
    tstart     float null,
    fids       varchar(20) null, 
    ascdsver   varchar(32) null, 
    ra     	 double precision null,
    dec    	 double precision null,
    roll   	 double precision null,
    process_status varchar(128) null,
    category_id int null,
    PRIMARY KEY (obsid)
);

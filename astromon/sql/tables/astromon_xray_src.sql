CREATE TABLE astromon_xray_src (
    obsid      int not null,
    id         int not null,
    name       varchar(22), -- ala 'CXOU J123456.7+765432'
    ra         double precision, 
    dec        double precision,
    net_counts float,  
    y_angle   float,  
    z_angle   float,  
    r_angle   float,  
    snr        float,
    near_neighbor_dist float,
    double_id int,
    status_id  int null,
    pileup float,
    acis_streak int,
    PRIMARY KEY (obsid,id)
);
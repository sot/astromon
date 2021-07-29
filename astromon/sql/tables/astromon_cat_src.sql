CREATE TABLE astromon_cat_src (
    obsid      int not null,
    id         int not null,
    catalog    varchar(16), 
    name       varchar(24),
    ra         double precision, 
    dec        double precision,
    mag        float null,
    y_angle   float,  
    z_angle   float,  
    PRIMARY KEY (obsid,id)
);
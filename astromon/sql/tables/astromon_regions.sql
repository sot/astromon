CREATE TABLE astromon_regions (
    region_id integer primary key autoincrement not null,
    ra    double precision not null,
    dec    double precision not null,
    radius    float not null,
    obsid int  not null,
    user varchar(50) not null,
    comments varchar(200)
);
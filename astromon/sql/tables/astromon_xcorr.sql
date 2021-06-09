CREATE TABLE astromon_xcorr (
    select_name varchar(14) not null,
    obsid int  not null,
    c_id  int not null,
    x_id  int not null,
    dy    float not null,
    dz    float not null,
    dr    float not null,
	PRIMARY KEY (select_name,obsid,c_id,x_id)
);
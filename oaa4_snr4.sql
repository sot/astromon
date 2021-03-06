-- Standard cross X-ray -- Catalog cross correlation query
--
--  X-ray sources within 3 arcmin off-axis (NONE) or 0.4 arcmin (grating)
--  X-ray snr >  5.0
--  X-ray extr_rad_grating,r,a,0.4,,,"Extr. rad. around grating source (arcmin)"
--
--  X-ray - Catalog position match within 3 arcsec
-- 
--  Catalog is relatively high precision:
--   'Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON', SDSS, 2MASS, USNO-B1.0
--
SELECT *, x.id  AS x_id,
          x.ra  AS x_ra,
          x.dec AS x_dec,
	  x.name AS x_name,
          c.id  AS c_id, 
          c.ra  AS c_ra,
          c.dec AS c_dec,
	  c.name AS c_name,
          (c.y_angle - x.y_angle) AS dy,
          (c.z_angle - x.z_angle) AS dz,
          (c.y_angle-x.y_angle)*(c.y_angle-x.y_angle) + (c.z_angle-x.z_angle)*(c.z_angle-x.z_angle) AS dr2
  FROM  astromon_xray_src AS x 
--  WHERE x.obsid NOT IN (SELECT unique obsid FROM astromon_xcorr)
    JOIN astromon_cat_src AS c ON x.obsid = c.obsid
    JOIN astromon_obs     AS o ON x.obsid = o.obsid
      WHERE o.process_status = 'OK'
	AND x.snr > 4.0
---     AND (x.status_id = NULL or x.status_id = 0) [KEEP so they can be seen if desired]
        AND x.near_neighbor_dist > 6.0
        AND (x.r_angle < 24 OR (o.grating = 'NONE' AND x.r_angle < 240))
        AND c.catalog IN ('Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON', 'SDSS', '2MASS', 'USNO-B1.0')
        AND (c.y_angle-x.y_angle)*(c.y_angle-x.y_angle) + (c.z_angle-x.z_angle)*(c.z_angle-x.z_angle) < 9.0


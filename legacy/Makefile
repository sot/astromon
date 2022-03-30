# Define the task name
TASK = astromon

SHARE = get_cat_obs_data.pl plot_offsets.pl xcorr.pl install_plots.pl
DATA = ASTROMON_table.rdb ICRS_tables astromon.par astromon_table_defs standard_xcorr.sql oaa4_snr4.sql task_schedule.cfg

# New-style location for data within share
INSTALL = $(SKA)
INSTALL_SHARE = $(INSTALL)/share/$(TASK)
INSTALL_DATA = $(INSTALL_SHARE)/data

test: 
	./get_cat_obs_data.pl -force 6893 8538

install:
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/

	mkdir -p $(INSTALL_DATA)/standard_xcorr
	if [ ! -e $(INSTALL_DATA)/standard_xcorr/www ] ; then \
           ln -s $(SKA)/www/ASPECT_PUBLIC/celmon $(INSTALL_DATA)/standard_xcorr/www ; fi

	mkdir -p $(INSTALL_DATA)/oaa4_snr4
	if [ ! -e $(INSTALL_DATA)/oaa4_snr4/www ] ; then \
           ln -s $(SKA)/www/ASPECT/celmon $(INSTALL_DATA)/oaa4_snr4/www ; fi

	mkdir -p $(INSTALL_SHARE)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/



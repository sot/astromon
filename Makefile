# Define the task name
TASK = astromon

# Define the installed executables for the task.  This directory is reserved
# for documented tools and not for dedicated scripts, one-off codes etc
# BIN = 

# Astromon uses this dedicated script to calculate offsets, accounting for
# improvements in boresight calibration since the data products were processed
SHARE = get_cat_obs_data.pl

# These are input data file required by astromon
DATA = ASTROMON_table.rdb ICRS_tables astromon.par astromon_table_defs

# Files that need to be installed in the local area for test purposes.
# Rules defined in Makefile.FLIGHT attempt first to locate files in the
# local t/ directory, followed by the flight root (/proj/sot/{ska,tst}). 
# In this case, the first two come from t/ while the scat program is
# from /proj/sot/ska
TEST_DEPS = bin/sun4/scat

# Set Flight environment to be SKA.  The other choice is TST.  Include the
# Makefile.FLIGHT make file that does most of the hard work
FLIGHT = SKA
include /proj/sot/ska/include/Makefile.FLIGHT

# Test target.  It is critical to have the "check_install" dependency to ensure
# that the install does not accidentally stomp on the flight version.  If this 
# does happen use /proj/sot/.snapshot to restore flight code immediately
#
# For this code, need to 'setenv SKA $PWD' for testing, since it requires an
# absolute path for SKA to function.  This should be fixed!

#	if [ -d $(INSTALL_DATA)/Obs_data/obs139 ] ; then rm -r $(INSTALL_DATA)/Obs_data/obs139 ; fi

test: check_install $(BIN) $(TEST_DEPS) install
	mkdir -p $(INSTALL_DATA)/Obs_data
	$(INSTALL_SHARE)/get_obs_xcorr_data.pl 6451

test_quick: check_install $(BIN) $(TEST_DEPS) install
	mkdir -p $(INSTALL_DATA)/Obs_data
	$(INSTALL_SHARE)/get_obs_xcorr_data.pl -force 139

install:
ifdef BIN
	mkdir -p $(INSTALL_BIN)
	rsync --times --cvs-exclude $(BIN) $(INSTALL_BIN)/
endif
ifdef DATA
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
endif
ifdef SHARE
	mkdir -p $(INSTALL_SHARE)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
endif


#!/usr/bin/env perl

use warnings;
use strict;

use IO::All;
use Getopt::Long;
use Ska::DatabaseUtil qw(:all);
use DBI;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use File::Basename qw(dirname);
use Cwd qw(abs_path);

our $ASTROMON_DATA  = abs_path(dirname(__FILE__)) . "/data";

our %opt = (select => 'standard_xcorr',
           );

GetOptions(\%opt,
           'select=s',
           'loud',
           )
  or die "Couldn't get options, stopped ";

my $xcorr_file;
for ("$opt{select}.sql", "$ASTROMON_DATA/$opt{select}.sql") {
    if (-r $_) {
        $xcorr_file = $_;
        last;
    }
}
die "Couldn't find $opt{select} file, stopped " unless $xcorr_file;

my %table_def = eval scalar io("$ASTROMON_DATA/astromon_table_defs")->slurp;

my $xcorr_query < io($xcorr_file);

# Connect to aca database with read/write access
my $dbh = DBI->connect("DBI:SQLite:dbname=$ASTROMON_DATA/astromon.db3", "", "");
die "Couldn't connect: $DBI::errstr" unless ($dbh);

# Do the actual query to cross-correlate the X-ray and Catalog data for each obsid
my @xcorr_data = sql_fetchall_array_of_hashref($dbh, $xcorr_query);

my @obsids = uniq map { $_->{obsid} } @xcorr_data;

# Delete any previous xcorr data for this select name
sql_do($dbh,
       qq/DELETE FROM astromon_xcorr WHERE select_name='$opt{select}'/);

# Insert the new data
my @astromon_xcorr_cols = keys %{ {grep {not ref} @{$table_def{astromon_xcorr}}} };

foreach my $xcorr_data (@xcorr_data) {
    printf("Inserting data for obsid=%d c_id=%d x_id=%d\n",
           $xcorr_data->{obsid},
           $xcorr_data->{c_id},
           $xcorr_data->{x_id},
          ) if $opt{loud};
    $xcorr_data->{select_name} = $opt{select};
    $xcorr_data->{dr} = sqrt($xcorr_data->{dr2});
    sql_insert_hash($dbh,
                    'astromon_xcorr',
                    { map { $_ => $xcorr_data->{$_} } @astromon_xcorr_cols }
                   );
}

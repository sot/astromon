#!/usr/bin/env python

import argparse

import Ska.DBI

parser = argparse.ArgumentParser(description="Convert sybase table to sqlite")
parser.add_argument("--table", default="timelines", help="Table name")
parser.add_argument("--sqlite-file", default="test.db3", help="Sqlite database file")
parser.add_argument("--schema-file", help="Schema file (default=auto-generated)")
parser.add_argument("--not-null", default=[], action="append", help="NOT NULL column")
parser.add_argument("--constraint", action="append", default=[], help="NOT NULL column")
args = parser.parse_args()


dbsy = Ska.DBI.DBI(dbi="sybase", server="sybase", user="aca_read")
dbs3 = Ska.DBI.DBI(dbi="sqlite", server=args.sqlite_file)


# Get table schema
cursy = dbsy.conn.cursor()
cursy.execute("SET CHAINED OFF")
cursy.callproc("sp_cols", [args.table])
cols = cursy.fetchall()

# Drop an existing table
try:
    dbs3.execute("drop table {}".format(args.table))
except Exception as err:
    print("Got an error dropping table: {}".format(err))

if args.schema_file:
    create = open(args.schema_file, "r").read()
else:
    constraints = ["CONSTRAINT {}".format(x) for x in args.constraint]
    not_nulls = set(x.lower() for x in args.not_null)
    fields = []
    for colname, coltype, collen in cols:
        colnull = " NOT NULL" if colname.lower() in not_nulls else ""
        fields.append("{} {}{}".format(colname, coltype, colnull))
    create = "CREATE TABLE {} (\n{})".format(
        args.table, ",\n".join(fields + constraints)
    )

print("Executing {}".format(create))
dbs3.execute(create)

print("Reading rows from sybase")
cursy.execute("SELECT * FROM {}".format(args.table))
rows = cursy.fetchall()
cursy.close()

print("Inserting {} rows".format(len(rows)))
curs3 = dbs3.conn.cursor()
insert = "insert into {} values ({})".format(args.table, ",".join(["?"] * len(cols)))
for i, row in enumerate(rows):
    curs3.execute(insert, row)
    print(i, "\r")

print("\n")
dbs3.commit()
curs3.close()

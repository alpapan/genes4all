#!/bin/bash
#$Id$

# Example File to start the daemon on your server
# For security, ensure that permissions for this file are chmod go-rwx and NOT owned by the Apache/Web user

# FULL PATH to Directory where biosoftware_bench files are produced (your Drupal 
# dir can have symlink pointing to this directory)

# Use current directory:
DIR_BENCH=$PWD
#or not (no / at the end):
DIR_BENCH="/data_storage/www-cluster/bench"

# User to run the daemon as; here is as current user. If it is another user, then you have to
# run this program as root (or sudo with su priviliges)
# NB: The user must have write privilages on DIR_BENCH
USER_BENCH=alexie
# Full path to daemon
DAEMON_BENCH=$DIR_BENCH/biosoftware_bench_daemon.pl
#DSN to connect to your Drupal database; see postgres config pg_hba.ini for setting password-less access but still provide a dummy password
DSN_BENCH='dbi:Pg:dbname=drupal_sra;host=speed;port=5433;user=www_db;password=3oioaWvsas1'

# go to dir and notify user we are there
cd $DIR_BENCH
if [ $PWD != $DIR_BENCH ]
  then
    echo "Directory $DIR_BENCH not found"
    exit
fi
echo "Will delete all biosoftware_bench output files in this directory, otherwise press control-C now:"
pwd
sleep 5

# stop previous server
$DAEMON_BENCH -dir $DIR_BENCH -user $USER_BENCH stop -dsn $DSN_BENCH
sleep 1

# delete all old files
rm -f *query *wait *done *error *log *submission *png *txt *html *gff *output *temp *xml *raw

# start new server. notice -cpus variable which you can control (used for BLASTALL)
$DAEMON_BENCH -dir $DIR_BENCH -user $USER_BENCH -cpus 3 -dsn $DSN_BENCH

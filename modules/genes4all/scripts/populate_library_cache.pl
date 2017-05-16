#!/usr/bin/perl -w

#unaligned query:
#select uniquename from feature where type_id=(select cvterm_id from cvterm where name='contig') and feature_id NOT IN (select feature_id from featureprop where type_id IN (select cvterm_id from cvterm where cv_id=(select cv_id from cv where name='ic_data_restrictions'))) and organism_id NOT IN (select organism_id from organismprop where type_id IN (select cvterm_id from cvterm where cv_id=(select cv_id from cv where name='ic_data_restrictions')));
# or
#select uniquename from feature where type_id=(select cvterm_id from cvterm where name='contig');
use strict;
$|=1;
my $cmd =`which wget`;chomp($cmd);
my $url = 'http://localhost';
my $csv = shift;

unless ($cmd && $url && $csv && -s $csv){die;}

$url .='/genes4all/info/feature?feature_id=';
$cmd .= ' -q -O /dev/null ';

open (IN,$csv);
my $header=<IN>;
chomp($header);
while (my $ln=<IN>){
	chomp($ln);
	last if $ln=~/^\(/;
	system($cmd.'"'.$url.$ln.'"');
	print '.';
	sleep(1); # to allow for stopping it
}

#!/usr/bin/perl -w

use strict;
my $cwd=`dirname $0`;chomp($cwd);
chdir($cwd);
my $pid;
open (PID,"blast_daemon.lock")||system("./start_daemon.sh");
while (my $ln=<PID>){
	if ($ln=~/^(\d+)/){
		$pid=$1;
		my $exists = kill 0, $pid;	#root or user who own daemon only
		chomp($exists);
		unless ($exists){
			warn "PID $pid was not found. Starting daemon anew\n";
			system("./start_daemon.sh");
		}
	}
}
close PID;
if (!$pid){
	warn "No lock file found $cwd/blast_daemon.lock Starting daemon anew\n";
        system("./start_daemon.sh");


}

;$Id:$

Drupal Bioinformatics Software Bench (biosoftware_bench)
============================================

To install

  * Install the Tabs module (the development version or a post-May 2010 stable version (if it exists) http://drupal.org/project/tabs).
  * Get Jquery Update 2-0 (currently an alpha release) from http://drupal.org/project/jquery_update
  * Get Jquery UI (stable version) from http://drupal.org/project/jquery_ui
  * Unpack above modules in your sites/all/modules/ directory
  * Read README.txt of Jquery_ui and download Jquery 1.7 as per instructions (1.6 is not available)
  * Go to Administer -> Site building -> Modules and enable the above modules
  ** Tabs may need the following settings to be set at /admin/settings/tabs due to tabs' bugs (at least pre-May 2010 stable version)
  ** navigation buttons: enabled ; descriptive tab URLs: enabled.
  * Install the gmod_dbsf module (stable version; http://drupal.org/project/gmod_dbsf)
  * place the biosoftware_bench directory into your modules directory.
  * Go to Administer -> Site building -> Modules and enable the module
  * Now go to Administer -> Site Configuration -> biosoftware_bench settings 
  * setup the BLAST server by first activating it and the formatdb/fastacmd (under "Setup the software available")
  * upload some databases. Some test databases which work with the demo are included under examples/ in this directory.

Daemon Script
======
A script to be used as a daemon is included. Disclaimer under GPL: Use them (and this module) at your own risk but note
that I 've been using them with http://www.insectacentral.org with no problems.

The daemon can fork with a specific unpriviliged user by specifying the -user argument. This function will require that you
run the script with root priviliges (the script will then change user and forks out). We use it at insectacentral.org like so
(some variables changed for security reasons):

biosoftware_bench_daemon.pl -dir /dbsf -user dbsf-user stop -dsn 'dbi:Pg:dbname=drupal-db;host=localhost;port=5432;user=www-db;password=123'
biosoftware_bench_daemon.pl -dir /dbsf -user dbsf-user -cpus 3 -dsn 'dbi:Pg:dbname=drupal-db;host=localhost;port=5432;user=www-db;password=123'

The password needs not be 'correct' if you have set connections to trust in your pg_hba file. The first command stops any running
server on that directory, and the second command start one. It will use 3 CPUs. The CPU argument is not used if you are using a 
Condor installation. The daemon checks if Condor is active at startup. 
So if you do activate Condor (via the Administer menu) then you need to restart the daemon.

How to autostart daemon:
------------------------

Two things need to be done to autostart daemon after system restart/shutdown. They are:
1. Editing start_daemon.sh
2. Adding check_daemon.pl script to crontab.

check_daemon.pl
---------------

check_daemon.pl script checks if the daemon is running. If not running, it starts a new one. This script makes use of start_daemon.sh. Find more details below.

start_daemon.sh
---------------

start_daemon.sh script stores config information to restart the daemon. You need to edit this script to set different variables. The following are the variables.

DIR_BENCH - temporary directory to which this program can write.
USER_BENCH - user name in whose account this program is going to run e.g. www-data.
DAEMON_BENCH - full address of biosoftware_bench_daemon.pl
DSN_BENCH - connection information

crontab
-------

sudo crontab -e

The above command can be used to edit crontab. To run the script every 5 minutes, add the following line.

5 * * * * /path/to/your/modules/biosoftware_bench/scripts/check_daemon.pl

Demonstration
=============
The files needed to run the demonstration link of BLAST are too large to distribute here. Please get them from
http://gmod-dbsf.googlecode.com/files/dbsf_demodbs.tar.bz2

The *.fsa (FASTA) files are sequence files used as input (query) sequences for the demonstration button of the DBSF BLAST server. We provide the
arabidopsis sequence files used as the subject databases. Feel free to use the .[pn]al or the formatted databases directly, but you will need to edit
the .[pn]al files first to make sure they have the full and correct path to the formatdb database.

 Then add the datasets at Administer -> Site Configuration
 -> Drupal Bioinformatic Server framework settings -> Setup available reference datasets/databases

Then set up the correct path for the BLAST databases under "Dataset variable". You must also provide at least one Group

Then add them as a dataset, the Demonstration expects the nucleotide file will use "Arabidopsis thaliana genes" 
as the friendly name and that the protein file will use "Arabidopsis thaliana proteins".


Regarding other modules
=============================
*  This module make use of the Tabs module. You may wish to add this to your tabs CSS to highlight selected tabs.
	.ui-tabs-nav li.ui-tabs-selected a,.ui-tabs-nav li.ui-tabs-active a, .ui-tabs-nav li.ui-tabs-selected a:hover, .ui-tabs-nav li.ui-tabs-active a:hover {
	background: #8ab9ff;
	color: #000000;
	}

Known bugs & incompatibilities
==============================
 * Tabs/CSS: internet explorer 7 has a bug with double row tabs. use IE8
 * Note that we are using the BLASTALL utility, not the new NCBI C++ Toolkit.
 
Maintainers
-----------
 Alexie Papanicolaou (alpapan)
 Temi Varghese (temivarg)

Other
=====
 If you find any bugs please report them to the issue tracker via http://drupal.org/project/biosoftware_bench

TODO list
 * Install and use the advanced_help module (http://drupal.org/project/advanced_help)

;$Id:$

GMOD Drupal Bioinformatic Server Framework (gmod_dbsf)
======================================================

This module is the first Drupal implementation for BioInformatics. It is GMOD (http://www.gmod.org) compatible
and aims to integrate with GMOD's Chado database schema (which is PostGres specific). For that reason, it requires
PostGres (tested with 8.3+) and make ample usage of foreign keys to ensure data integrity. A number of tables prefixed
with gmod_dbsf_ will be created in your Drupal database and they will follow the Chado conventions.

Disclaimer:
A number of perl scripts are found in the scripts/ sub-directory. Use them (and this module) at your own risk but note 
that I 've been using them with http://www.insectacentral.org with no problems. All scripts have perldocs.

To install
==========
  
  * You must be using a PostGres database engine to use this module (currently) 
  * Make sure BioPerl is installed.
  * Install all Perl modules needed (see below)
  * Optionally, install Chado (http://www.gmod.org) including the gmod_* scripts. Currently this will do the trick
    svn co https://gmod.svn.sourceforge.net/svnroot/gmod/schema/trunk # and check out chado/
  * We included an improved version of gmod_bulk_load_gff3 to allow for less cumbersome organism databasing
  * Optionally, Make a Chado database and edit your Drupal site settings.php (e.g. drupal/sites/default/settings.php) according to advice in next section
  * Recommended: Try running every Perl script in the scripts sub-directory to ensure they don't produce an error for your installation.
    Some scripts will expect to use Chado, don't use them if you are not installing Chado
  * Unless already present in your distribution (includes/dynatree/jquery.dynatree.min.js), download the dynatree javascript library; see includes/dynatree/README
  * place the gmod_dbsf folder into your modules directory. Either in sites/all/modules (available then globablly) or a specific site 
  * Go to Administer -> Site building -> Modules
    ** Enable the module under gmod_dbsf 
    ** Note, I recommend you enable the modules sequentially according to dependencies. I have experienced problems when i tried to install multiple modules
    at the same time.
  * This will install a number of database tables prefixed with genes4all_

To Uninstall
============

 *** Please first disable all dependent modules
 *** Uninstall /them/ 
 *** Then de-activate and unistall this base module.


settings.php
============
This is the file which sets all the basic settings from Drupal to work. First we recommend the following options if they don't exist already

ini_set('memory_limit', '128M');
ini_set('upload_max_filesize','100M');
ini_set('post_max_size','100M');

Also set cookie and base_url, for example:
$base_url = 'http://insectacentral.org';  // NO trailing slash!
$cookie_domain = 'insectacentral.org';

If you are using the chado module, find and edit the $db_url array. It defines which databases are accessible by Drupal

$db_url = array (
        'default'=>'pgsql://<user with read/write access>:<password>@<external server hostname or localhost for local>:<db port>/<name of drupal db>',
// this is needed only for chado. this username has only read access and is used for most operations
       'chado'=>'pgsql://<user with read access to chado>:<password>@<external server hostname or localhost for local>:<db port>/<name of chado db>',
// we use a different username for write access, that is more secure:
        'chado_edit'=>'pgsql://<user with read/write access to chado>:<password>@<external server hostname or localhost for local>:<db port>/<name of chado db>',
//these are only for genes4all_explore which is a separate sub-module of the genes4all module
       'go'=>'pgsql://<user with read access>:<password>@<external server hostname or localhost for local>:<db port>/<go db if the genes4all_explore is present>',
        '77259'=>'pgsql://<user with read access>:<password>@<external server hostname or localhost for local>:<db port>/<Bio::SeqFeature::Store db for species with NCBI TaxID 77259>',
        );

You have to replace the <values> with the correct credentials for your setup. Make sure that no other $db_url variable is declared in your settings.php file

Note that the chado databases should NOT have a prefix; e.g.
$db_prefix = array (
        'default' =>'',
        'chado'=>''
        );
is OK

$db_prefix = array (
        'default' =>'mydrupal_',
        'chado'=>''
        );
is OK

$db_prefix = array (
        'default' =>'mydrup_',
        'chado'=>'chado'
        );

is NOT OK

Sub-modules
===========

If any submodules exist in the modules/ directory, please read the README in their subdirectory as well.
They are activated separately but you don't need to move them (keep them as a submodule).

Perl modules needed
===================
Getopt::Long
Pod::Usage
Time::Progress
DBI
BioPerl (http://www.bioperl.org)
Bio::Graphics;


Maintainers
-----------
Alexie Papanicolaou (alpapan)

Other & bugs
============
 
 If you find any bugs please report them to the issue tracker via http://drupal.org/project/gmod_dbsf

TODO list
==========
 * Install and use the advanced_help module (http://drupal.org/project/advanced_help)

;$Id:$

Genes4all, visualization modules for -omics of non-model species
================================================================

To install

  * Download the Dojo Toolkit (tested with 1.5) from http://dojotoolkit.org/download/ Make it is the 'Dojo Toolkit Release' version.
  * Unpack it in this directory and, if needed, rename the directory dojo-* to lib so that lib/dijit, lib/dojo and lib/dojox are present.
  ** (FAQ: Do I need the entire dojo Toolkit? No but until we stabilize which elements are needed, please just unpack all of it)
  * Install the Tabs module (the development version or a post-May 2010 stable version (if it exists) http://drupal.org/project/tabs).
  * Get Jquery Update 2-0 (currently an alpha release) from http://drupal.org/project/jquery_update
  * Get Jquery UI (stable version) from http://drupal.org/project/jquery_ui
  * Unpack above modules in your sites/all/modules/ directory
  * Read README.txt of Jquery_ui and download Jquery 1.7 as per instructions (1.6 is not available)
  * Go to Administer -> Site building -> Modules and enable the above modules
  ** Tabs may need the following settings to be set at /admin/settings/tabs due to tabs' bugs (at least pre-May 2010 stable version)
  ** navigation buttons: enabled ; descriptive tab URLs: enabled.
  * Install the gmod_dbsf module (stable version; http://drupal.org/project/gmod_dbsf)
  * place the genes4all directory into your modules directory.
  * Go to Administer -> Site building -> Modules and enable the module

Regarding other modules
=============================
*  This module makes use of the Tabs module. You may wish to add this to your tabs CSS to highlight selected tabs.
	.ui-tabs-nav li.ui-tabs-selected a,.ui-tabs-nav li.ui-tabs-active a, .ui-tabs-nav li.ui-tabs-selected a:hover, .ui-tabs-nav li.ui-tabs-active a:hover {
	background: #8ab9ff;
	color: #000000;
	}
* genes4all_curate: will need the biosoftware_bench project
* genes4all_curate: See README in module/genes4all_curate/


Known bugs & incompatibilities
==============================
 * Tabs/CSS: internet explorer 7 has a bug with double row tabs. use IE8
 * biosoftware_bench: Note that we are using the BLASTALL utility, not the new NCBI C++ Toolkit.
 
Maintainers
-----------
 Alexie Papanicolaou (alpapan)

Other
=====
 If you find any bugs please report them to the issue tracker via http://drupal.org/project/genes4all

TODO list
 * Install and use the advanced_help module (http://drupal.org/project/advanced_help)

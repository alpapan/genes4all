;$Id:$

JBrowse module for Drupal
==========================

This module requires a re-written Chado perl adaptor: Bio::DB::Das::Chado_AP and Bio::DB::Das::Chado_AP::Segment.pm
It can be installed site wide (in perl library dir), installed for the Apache user only
or (I think) provided in JBrowse directory: lib/Bio/DB/Das/Chado_AP.pm and lib/Bio/DB/Das/Chado_AP/*

Tested with JBrowse 1.2

This will module will probably not work on systems with SELINUX (unless it is configurated appropriately).

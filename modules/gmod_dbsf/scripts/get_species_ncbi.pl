#!/usr/bin/perl -w
#$Id$

=pod

=head1 NAME

  get_species_ncbi.pl

=head1 USAGE

 'species:s' => Species or common name
 'genus:s'  => Genus
 'ncbi:i'   => NCBI Taxonomy ID
 'all'	    => Return entire phylogeny
 'd|flatfile_dir:s'  => Bio::DB::Taxonomy flatfile (optional)

=head1 DESCRIPTION

 An alternative is to use a flatfile acquired from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=head1 AUTHORS

    Alexie Papanicolaou

    Centre for Ecology and Conservation, University of Exeter, UK
    * recipient of complaints:  alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

 This software is released under the GNU General Public License version 3 (GPLv3).
 It is provided "as is" without warranty of any kind.
 You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
 Please note that incorporating the whole software or parts of its code in proprietary software 
 is prohibited under the current license.


=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::Taxon;
use Bio::DB::Taxonomy;

my ($genus,$species,$ncbi_taxid,$all,$verbose);
my $flatfile_dir='/db/ncbi_taxonomy/';
GetOptions(
 'species:s' => \$species,
 'genus:s'  => \$genus,
 'ncbi:i'   => \$ncbi_taxid,
 'all'	   => \$all,
 'd|flatfile_dir:s' => \$flatfile_dir,
 'verbose' => \$verbose,
);
if ($ncbi_taxid && $ncbi_taxid!~/^\d+$/){$ncbi_taxid='';}

unless ($ncbi_taxid || $genus || $species){
	$ncbi_taxid=shift;
	if ($ncbi_taxid!~/^\d+$/){$genus=$ncbi_taxid;$species=shift;$ncbi_taxid='';}
	warn ("Need an NCBI id or both a genus and species\n") unless ($ncbi_taxid || ($genus && $species));
	pod2usage unless ($ncbi_taxid || ($genus && $species));
}
my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
   warn "Using $flatfile_dir as NCBI TaxDB directory\n" if $verbose;
   $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile',-nodesfile => $flatfile_dir.'/nodes.dmp',-namesfile => $flatfile_dir.'/names.dmp',-directory => $flatfile_dir);
}
if ($genus && $species && !$ncbi_taxid){
    $ncbi_taxid = $taxondb->get_taxonid($genus.' '.$species);
    if (!$ncbi_taxid){
        die "Could not find NCBI taxid for $genus $species.\n";
    }
}
my (%all_taxa,%all_taxa_sorted);
my $sorted_ref=\%all_taxa_sorted;
my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your query\n";;
# we want class, order, family, genus, species, common name
my ($class,$order,$family,$i);
 my $species_rank=$taxon->rank();
 $species = $taxon->scientific_name() if ($species_rank eq 'species');
 if (!$species){
	my $t=$taxon;
	while ($t=$t->ancestor()){
		$species_rank=$t->rank();
		next if !$species_rank;
		$species = $taxon->scientific_name() if ($species_rank eq 'species'); 
	}
	if (!$species){
	 	die "Cannot find a scientific species name for your query at $species_rank\n";
	}
 }
 $species =~s/^(\w+)\s//;
 $i++;$sorted_ref->{$i}={
	'rank' => 'species',
	'name' => $species,
 } if $all;
 $all_taxa{$species_rank}=$species if $all;
 my $t=$taxon;
 while ($t=$t->ancestor()){
        my $rank=$t->rank();
	next if !$rank;
	my $name = $t->scientific_name();
	next if !$name;
	$i++;$sorted_ref->{$i}={
		'rank' => $rank,
		'name' => $name,
	} if $all;
	$all_taxa{$rank}=$name;
	if ($rank eq 'genus'){
	 	$genus=$name;
		last;
	}
 }
if (!$genus || !$species){die "Cannot find genus or species\n";}
 my $common=$taxon->common_names() ? $taxon->common_names() : $genus.' '.$species;
 while ($t=$t->ancestor()){
 	my $rank=$t->rank();
	next if !$rank;
	my $name = $t->scientific_name() ;
	next if !$name;
	$i++;$sorted_ref->{$i}={
		'rank' => $rank,
		'name' => $name,
	} if $all;
	$all_taxa{$rank}=$name;
 	if ($rank eq 'class'){
 		$class=$name;
 		last if (!$all);
 	}elsif($rank eq 'order'){
 		$order=$name;
 	}elsif($rank eq 'family'){
        	$family=$name;
	}
 }
if (!$all){
	if (!$class){
		# go up
		if ($all_taxa{'phylum'}){
			$class=$all_taxa{'phylum'};
		}elsif ($all_taxa{'kingdom'}){
			$class=$all_taxa{'kingdom'};
		}else{
			$class='unknown';
		}
	}if (!$order){
		#go down one, then up one
		if ($all_taxa{'suborder'}){
			$order=$all_taxa{'suborder'};
		}elsif ($all_taxa{'subclass'}){
	                $order=$all_taxa{'subclass'};
	        }else{
	                $order='unknown';
	        }
	}if (!$family){
		#go down one, then up one
		if ($all_taxa{'subfamily'}){
	                $family=$all_taxa{'subfamily'};
        	}elsif ($all_taxa{'superfamily'}){
	                $family=$all_taxa{'superfamily'};
        	}else{
	                $family='unknown';
        	}
	}
	print "class:$class;order:$order;family:$family;genus:$genus;species:$species;common:$common;ncbi:$ncbi_taxid\n";
}else{
	my ($print);
	foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
		$print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';';
	}
	print $print."common:$common;ncbi:$ncbi_taxid\n";
}

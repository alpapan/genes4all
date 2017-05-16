#!/usr/bin/perl -w
#$Id$

=pod

=head1 NAME

 ic_gff2all.pl

=head1 USAGE

 'infile:s' => \$infile
 'genbank'=> \$genbank,
 'fasta' => \$fasta,

=head1 AUTHORS

 Alexie Papanicolaou 1 2

    1 Max Planck Institute for Chemical Ecology, Germany
    2 Centre for Ecology and Conservation, University of Exeter, UK
    alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind. You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. Please note that incorporating the whole software or parts of its code in proprietary software is prohibited under the current license.

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::Annotation::Collection;
use Data::Dumper;
use Bio::Taxon;
use Bio::DB::Taxonomy;


my ($infile,$genbank,$fasta,$species,$ncbi_taxid,$compress);
my ($fasta_out,$gb_out);
my $authority='insectacentral.org';
my $gzip;
my %species_hash;
GetOptions(
 'infile:s' => \$infile,
 'genbank'=> \$genbank,
 'fasta' => \$fasta,
 'species:s'=> \$species,
 'authority:s'=>\$authority,
 'compress' => \$compress,
);

unless ($infile && -s $infile){warn "Specify an input GFF3 file\n";pod2usage;}
unless ($genbank || $fasta){warn "Select at least one output format\n";pod2usage;}
if ($compress){
	$gzip=`which gzip`||die ("Cannot find gzip\n");
	chomp($gzip);
}
my $outfile=$infile;
$outfile=~s/\.gff3$//;
my $outfiles=$infile;
my $gff_in = Bio::Tools::GFF->new(-file => $infile, -gff_version => 3);

#TODO WARNING REQUIRED PATCHING /usr/local/share/perl/5.10.0/Bio/Tools/GFF.pm line 323
$gff_in ->features_attached_to_seqs(1);
if ($fasta){
    $fasta_out = Bio::SeqIO->new(-file => ">$outfile.fsa", format => 'fasta');
    $outfiles.=" $outfile.fsa";
}if ($genbank){
    $gb_out = Bio::SeqIO->new(-file => ">$outfile.gb", format => 'genbank');
    $outfiles.=" $outfile.gb";
}
my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
while(my $feature = $gff_in->next_feature()) {
	if ($feature->frame() eq "."){$feature->frame(0);}
	
}
my @seq_objs=$gff_in->get_seqs();
foreach my $seq_obj (@seq_objs){
    if ($genbank){
    	my $id=$seq_obj->id();
	    my ($data_type,$db_id,$assembly_v,$serial);
	    ## specific to est2assembly naming scheme
    	if ($id=~/^([A-Z]{2})(\d+)([A-Z][a-z])([A-Z][a-z]+)(\d+)/){
		$db_id=$1;
        $ncbi_taxid=$2;
		$assembly_v=$3;
		$data_type=$4;
		$serial=$5;
    	}
    	if ($ncbi_taxid && !$species){$species=$ncbi_taxid;}
	if ($data_type){
		if ($data_type eq "Apep"){
			$seq_obj->alphabet('protein');
		}else{
			$seq_obj->alphabet('dna');
		}		
		$seq_obj->authority($authority);
		$seq_obj->object_id($id);
		$seq_obj->namespace($db_id.$ncbi_taxid.$assembly_v.$data_type);
	}
	## END est2assembly naming scheme
        if ($species && !$ncbi_taxid){
            $ncbi_taxid = $taxondb->get_taxonid($species);
        }
        if (!$ncbi_taxid){
        	warn "Could not find NCBI taxid for $id.\n";
        }
        my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid) if ($ncbi_taxid && !$species_hash{$ncbi_taxid});
        $species_hash{$ncbi_taxid}=$taxon if ($ncbi_taxid && $taxon);
        $seq_obj->species($taxon);        
        #my $ann_col = new Bio::Annotation::Collection;$seq_obj->annotation($ann_col);
        #$seq_obj->get_SeqFeatures();
        $gb_out->write_seq($seq_obj);
    }

    if ($fasta){
        my $seq=$seq_obj->seq();
        $seq=~s/\-+//g;
        $seq_obj->seq($seq);
        $fasta_out->write_seq($seq_obj );
    }
}
$gff_in->close();
if ($gzip && $outfiles){
	print $outfiles."\n";
	system ("$gzip -9 $outfiles");
}

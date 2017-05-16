#!/usr/bin/perl -w
# $Id$

=pod

=head1 NAME

 index_fasta.pl
 
=head1 USAGE

 Index and/or retrieve a sequence from a FASTA flatfile using BioPerl
 
 To index a database, simple give -d <name of file>
 To retrieve, give -d <name of file> -s one or more ID/accessios
 
 * -d|database:s => name of databae
   -s|a|id:s        => ID accession to fetch
   -f|force            => force re-indexing

=cut
use strict;
use Bio::Index::Fasta;
use Text::Wrap;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
$Text::Wrap::columns = 80;
my $infile;
my @ids;
my ($force);
GetOptions(
			'd|infile:s'         => \$infile,
			'force'              => \$force,
			's|accession|id:s{,}' => \@ids,
);
pod2usage if !$infile;
die "Cannot find $infile" if !-s $infile;
my $inx;

if ($force) {
	unlink("$infile.index");
}
$inx = Bio::Index::Fasta->new( -filename => "$infile.index", -write_flag => 1 )
  if ( !-s "$infile.index" && -s $infile );
$inx->make_index($infile) if $inx;
$inx = Bio::Index::Fasta->new( -filename => "$infile.index" );
if (@ids) {
	my %hash; # to get them only once each
	foreach my $id ( sort @ids ) {
		$hash{$id}=1;
	}
	foreach my $id (keys %hash){
		my $seq_obj = $inx->fetch($id);
		if ($seq_obj) {
			my $seq = $seq_obj->seq();
			print ">" . $seq_obj->id() . "\n" . wrap( '', '', $seq ) . "\n"
			  if $seq;
		}else{
			print "\nFailed to find sequence $id\n";
		}
	}
}

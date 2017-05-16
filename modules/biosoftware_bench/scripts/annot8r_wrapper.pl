#!/usr/bin/perl -w
#$Id$

=head1 NAME

 annot8r_wrapper.pl - automating ANNOT8R for web-server

=head1 USAGE
 
 automating ANNOT8R for web-server

 
			'annot:s'     => \$annot8r_exec,
			'physprop:s'  => \$physprop_exec,
			'annot2gff:s' => \$bb_annot8r2gff_exec,
			'blast:s'     => \$blast_exec,
		*	'infile:s'    => \$infile,
			'seq:s'       => \$inseq,                 #allow to pass a string?
			'outfile:s'   => \$outfile,
			'dbgo:s'        => \$gofile,
			'dbec:s'        => \$ecfile,
			'dbkegg:s'      => \$keggfile,
		**	'go'      => \$dogo,
		**	'ec'      => \$doec,
		**	'kegg'      => \$dokegg,
			'blastout:s'  => \$blastout

=head1 AUTHORS

 Alexie Papanicolaou 1 2

        1 Max Planck Institute for Chemical Ecology, Germany
        2 Centre for Ecology and Conservation, University of Exeter, UK
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

None known so far. Please report them.


=cut

use strict;
use Pod::Usage;
use Getopt::Long;
$| = 1;

#default values for globals
my $physprop_exec       = "/home/alexie/bin/annot8r_physprop.pl";
my $annot8r2gff_exec = "/home/alexie/bin/ic_annot8r2gff.pl";
my $annot8r_exec        = "/home/alexie/bin/annot8rAP-gff.pl";
my $blast_exec          = "/usr/bin/blastall";
my $gofile              = "/db/blastdb/uniprot_GO.fsa";
my $ecfile              = "/db/blastdb/uniprot_EC.fsa";
my $keggfile            = "/db/blastdb/uniprot_KEGG.fsa";
my $blastout;    # if we already have a blastfile
my ($dogo,$doec,$dokegg);
my $inseq;       # a string
my $pgport=5436;
my $pghost='localhost';
my $cpus=2;
my ($infile,$query_count,$pguser,$pgpass,$blastmat,@gffinfiles,$organism,@gffoutfiles,$outfile);
GetOptions(
			'annot:s'     => \$annot8r_exec,
			'physprop:s'  => \$physprop_exec,
			'annot2gff:s' => \$annot8r2gff_exec,
			'blast:s'     => \$blast_exec,
			'matrixblast:s' => \$blastmat,
			'i|infile:s'    => \$infile,
			'seq:s'       => \$inseq,                 #allow to pass a string?
			'dbgo:s'        => \$gofile,
			'dbec:s'        => \$ecfile,
			'dbkegg:s'      => \$keggfile,
			'go'      => \$dogo,
			'ec'      => \$doec,
			'kegg'      => \$dokegg,
			'blastout:s'  => \$blastout,	# one search only!
			'host|pghost:s'   => \$pghost,
			'port|pgport:s'   => \$pgport,
			'user|pguser:s'   => \$pguser,
			'pass|pgpass:s'    => \$pgpass,
			'organism:s'     => \$organism,
			'o|output:s'     => \$outfile,
			'cpus:i'         => \$cpus,
);
if ( $infile && $inseq ) {
	die("Not both input file and raw\n"); 
}
if (!$infile && !$inseq){
	$infile =shift || warn ("Need an input file!\n");
}
unless (    
		 $annot8r_exec
		 && $physprop_exec
		 && $blast_exec
		 && ( ( $infile && -s $infile ) || $inseq ) )
{
	pod2usage;
	die;
}

$ENV{'BLASTMAT'}=defined $blastmat ? $blastmat : $ENV{'BLASTMAT'};
if (!$ENV{'BLASTMAT'}){
	die "No BLAST MATRIX directory given\n";
} 
my $dsn =" -pgport $pgport -pghost $pghost -pguser $pguser";
if ($pgpass){$dsn.=" -pgpass $pgpass ";}
#warn "Using $dsn in order to connect...\n";

#safety check
if ($infile ) { 
	($infile,$query_count) = &clean_file( $infile, "fasta", "bioperl" ); 
}
elsif ($inseq) {
	die("String passing not supported yet\n");
	$inseq = uri_escape($inseq);
	open( IN, ">user_input" );    # add some kind of job id?
	print IN $inseq;
	close IN;
	$infile = "user_input";
}

$outfile = defined $outfile ? $outfile : "$infile.annot8r.gff";
$dogo=1 if ($gofile && (-s $gofile.'.pin' || -s $gofile.'.pal'));
$doec=1 if ($ecfile && (-s $ecfile.'.pin') || -s $ecfile.'.pal');
$dokegg=1 if ($keggfile && (-s $keggfile.'.pin' || -s $keggfile.'.pal'));
#TODO: nucs -> blastx and get tile??
print "Running annot8r_wrapper for $query_count sequences\n";
print "Physical properties...\n";
#output is .phys
system("$physprop_exec $infile >/dev/null 2>$infile.error");

if ($dogo){
	print "Gene Ontology:";
	my $local_blast;
	if ($blastout) {$local_blast=$blastout; }
	else {$local_blast = "$infile.GO.blast"; }
	if ( -s $local_blast) {
		print "BLAST output $local_blast already exists. skipping blast.\n";}
	else {
		print " BLAST";
		system("$blast_exec -a $cpus -p blastp -i $infile -d $gofile -e 1e-20 -b 30 -v 30 -o $local_blast  2>>$infile.error >/dev/null"		);
		unless ( -s $local_blast) { warn("Cannot find GO blast output\n"); }
	}
	print " annot8r\n";
	system("$annot8r_exec -r GO $local_blast -skip -o $infile.annot8r.GO $dsn >/dev/null  2>>$infile.error");
	if ( !-s "$infile.annot8r.GO" ) { warn("Cannot find GO annot8r output\n"); }
	else{unlink($local_blast);}
	push(@gffinfiles,"$infile.annot8r.GO");
}
if ($doec){
    print "Enzyme classification:";
	my $local_blast;
	if ($blastout) {$local_blast=$blastout; }
	else {$local_blast = "$infile.EC.blast"; }
	if ( -s $local_blast) {print "BLAST output $local_blast already exists. skipping blast.\n";}
	else {
		print " BLAST";
		system("$blast_exec -a $cpus -p blastp -i $infile -d $ecfile -e 1e-20 -b 30 -v 30 -o $local_blast >/dev/null  2>>$infile.error"		);
		unless ( -s $local_blast) { warn("Cannot EC find blast output\n"); }
	}
	print " annot8r\n";
	system("$annot8r_exec -r EC $local_blast -skip -o $infile.annot8r.EC $dsn >/dev/null  2>>$infile.error");
	if( !-s "$infile.annot8r.EC") { warn("Cannot find EC annot8r output\n"); }
	else{unlink($local_blast);}
	push(@gffinfiles,"$infile.annot8r.EC");
}
if ($dokegg){
	print "KEGG terms:";
	my $local_blast;
	if ($blastout) {$local_blast=$blastout; }
	else {$local_blast = "$infile.KEGG.blast"; }
	if ( -s $local_blast) {print "BLAST output $local_blast already exists. skipping blast.\n";}
	else {
		print " BLAST";
		system("$blast_exec -a $cpus -p blastp -i $infile -d $keggfile -e 1e-20 -b 30 -v 30 -o $local_blast >/dev/null  2>>$infile.error");
		unless ( -s $local_blast) { warn("Cannot KEGG find blast output\n"); }
	}
	print " annot8r\n";
	system("$annot8r_exec -r KEGG $local_blast -skip -o $infile.annot8r.KEGG $dsn >/dev/null  2>>$infile.error");
	if ( !-s "$infile.annot8r.KEGG") { warn("Cannot find KEGG annot8r output\n"); }
	else{unlink($local_blast);}
	push(@gffinfiles,"$infile.annot8r.KEGG");
}
if (@gffinfiles){
	print "Preparing GFF files...\n";
	if ($organism){
		system("$annot8r2gff_exec -org $organism -nodb -p $infile.phys @gffinfiles 2>/dev/null >/dev/null");
	}else{
        system("$annot8r2gff_exec -nodb -p $infile.phys @gffinfiles  2>/dev/null >/dev/null");
	}
}

push (@gffoutfiles,"$infile.phys.gff") if -s "$infile.phys.gff";
unlink("$infile.phys") ;
foreach my $in (@gffinfiles){
	push (@gffoutfiles,"$in.gff") if -s "$in.gff";
	unlink $in;
}

open (OUT, ">$outfile");
print OUT "##gff-version 3\n";
foreach my $gff (@gffoutfiles){
	open (GFF,$gff);
	while (my $line=<GFF>){
		print OUT $line if ($line!~/^#/);
	}
	close (GFF);
	unlink $gff;
}
close (OUT);
print "Done\n";
############################################
sub clean_file($$) {
        my $in     = shift || die;
        my $format = shift || die;
        my $method = shift || die;
        my $count=int(0);
        my $keep   = "^A-Za-z0-9\-_.!~*'()\\n";
        my $out    = uri_escape($in) . ".cleaned";
        if ( $method =~ /bioperl/i ) {
                use Bio::SeqIO;
                my $filein  = new Bio::SeqIO( -file => $in,     -format => "$format" );
                my $fileout = new Bio::SeqIO( -file => ">$out", -format => "fasta" );
                while ( my $seq = $filein->next_seq() ) {
                        $count++;
                        #strip description.
                        #$seq->description("");
                        $fileout->write_seq($seq);
                }
        } elsif ( $method =~ /uri/i ) {
                use URI::Escape;
                open( IN,  $in );
                open( OUT, ">$out" );
                while ( my $line = <IN> ) {
                        
                        $line = uri_escape( $line, $keep );
                        print OUT $line;
                }
                close IN;
                close OUT;
        }
        system ("mv -f $out $in");
        return ($in,$count);
}

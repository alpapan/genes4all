#!/usr/bin/perl -w
# $Id$

=pod

=head1 NAME

  biosoftware_bench_daemon.pl

=head1 USAGE

  ** are required, * are recommended

     'help'       => Show usage options        
  ** 'dsn|database_connection:s'     => DSN connection file or text string for Drupal DB
  *  'dir:s'      => Directory to look for/create files. If condor, then it must be on a SHARED directory (def. /tmp)
  *  'user:s'     => Run as user (e.g. apache user - e.g. www-data). Needs root priviliges
     'dbprefix:s' => Drupal database prefix (see settings.ini)
     'cpus:i'     => \Number of CPUS if no condor
     'debug'    => Debug,
     'link|db_link:s'  => DB link for HMTL blast reports
     'stop'     => Kill existing daemon running in dir
     'status'   => Report if a daemon is already running in dir
     
 
	You can check status or stop existing daemons by giving 'status' or 'stop', respectively
	If you're root, you can set effective user with -user (defaults to existing user). 
	If using it for biosoftware_bench, we recommend setting it  to your apache user.
	Note that all options can abbreviated with the first unique set of letters

=head1 TIPS

 a dsn file is simple text file which has one line like this:
 dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321

 Currently, only one environmental variable per software is allowed.


=head1 CAUTION

 If using InterProScan, you must edit the converter.pl file with this change for HTML:
 
   my $xfile = $iprscan->getConfigValue('toolxml');
   die '"', __FILE__, '"', ", line \"", __LINE__, "\" No toolxml configured" unless $xfile;
 + system("$0 -format xml -input $input -job $jobid -out $xfile") if !-f $xfile;
   die '"', __FILE__, '"', ", line \"", __LINE__, "\" $xfile does not exist" unless(-f $xfile);
 

=head1 CHANGELOG

  01Feb10
    Added support for 
      * multiple lines in .wait file
      * iprscan (but not graphmaking)
   21Feb10
      * annot8r (but not graphmaking)

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
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use POSIX qw(setsid setuid);
use Fcntl ':flock';
use DBI;
use Bio::Tools::GFF;
use File::Basename;
use Bio::SearchIO;
use Bio::SearchIO::Writer::TextResultWriter;
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::SearchIO::Writer::GbrowseGFF;
use Bio::Graphics;
use Bio::FeatureIO;
use Bio::SeqFeature::Generic;
$| = 1;

# globals
#### my $condor_submit=`which condor_submit`;chomp($condor_submit);
my $dir        = '/tmp';
my $exist_user = `whoami`;
chomp($exist_user);

# so that we can revert if needed
my ( $exist_uid, $exist_gid ) = &get_uid($exist_user);
my $user = $exist_user;
my ( $status, $stop, $condor, $help, %software_data, $software_data_str, $logfile, $debug );
my $dbprefix   = ' ';
my $cpus       = 2;
my $dsn        = 'dbi:Pg:dbname=ic-drupal;host=localhost;port=5436;user=www-db;password=123';
my $db_link    = 'http://insectacentral.org/genes4all/search/accession?accession=%s';
my $sleep_time = 30;
GetOptions(
			'dsn:s'          => \$dsn,
			'help'           => \$help,
			'dir:s'          => \$dir,
			'user:s'         => \$user,
			'dbprefix:s'     => \$dbprefix,
			'cpus:i'         => \$cpus,        # if no condor
			'debug'          => \$debug,
			'link|db_link:s' => \$db_link,     # for HMTL blast reports
			'stop'           => \$stop,
			'status'         => \$status,
			'sleep:i'        => \$sleep_time
);
if ( !$dsn || !$user || !$dir || $help ) { pod2usage; }
my $lock = $dir . '/blast_daemon.lock';
my $action = shift if $ARGV[0];
## If capturing a signal, unlock
use sigtrap 'handler' => \&unlock, 'normal-signals';
umask 0;
if ($stop)   { &killme();   exit; }
if ($status) { &statusme(); exit; }
if ($action) {

	if ( $action =~ /stop/ ) {
		&killme();
		exit;
	} elsif ( $action =~ /status/ ) {
		&statusme();
		exit;
	} elsif ( $action =~ /start/ ) {
	} else {
		warn "Action can only be one of stop or status\n";
		pod2usage;
	}
}
print "\n";
opendir( DIR, $dir ) || die("Cannot open $dir\n");
my $escape_pattern = '[\$\^\&\*\(\)\<\\\?\:\+]+';
my $var_hash       = &get_pars();
my ( $new_uid, $new_gid ) = &change_uid($user);
&daemonize();
print "Will check for:\n";

if ( $var_hash->{'blastall'}->{'exec'} ) {
	print "\tblastall\n";
}
if ( $var_hash->{'ssaha2'}->{'exec'} ) {
	print "\tssaha2\n";
}
if ( $var_hash->{'ipr_convert'}->{'exec'} ) {
	print "\tInterPro\n";
}
if ( $var_hash->{'annot8r'}->{'exec'} ) {
	print "\tannot8r\n";
}
print "\n";

# Here is what the daemon actually does
while (1) {
	&check_submissions();
	sleep(5);    # I/O Friendship
	&check_results();
	sleep($sleep_time) unless $debug;
}

# before exiting while, unlock everything
&unlock();
###
sub daemonize() {
	chdir($dir) or die "Can't chdir to $dir: $!";
	if ( -s $lock ) {
		open( LOCK, $lock );
		while ( my $line = <LOCK> ) {
			if ( $line =~ /(\d+)/ ) {
				my $pid     = $1;
				my $failure = `ps -A |grep $pid`;
				if ( $failure =~ /^\s*(\d+)/ ) {
					die "Program already running as pid $1.\n";
				}
			}
		}
		close LOCK;
	}

	#two children down
	defined( my $pid = fork ) or die "Can't fork: $!";
	exit if $pid;
	defined( $pid = fork ) or die "Can't fork: $!";
	exit if $pid;
	setsid or die "Can't start a new session: $!";
	print "Started as pid $$\n";
	open( LOCK, ">$lock" ) or die "Cannot lock - $!\n";
	## pid of child (now master) -> write it in a lock file.
	print LOCK $$ . "\n";
	my $locked = flock( LOCK, LOCK_EX ) or die "Cannot lock - $!\n";
	$logfile = "$dir/$$.log";
	open( LOG, ">$logfile" ) || die("Cannot write in $dir");
	open STDIN,  '/dev/null'  or die "Can't read /dev/null: $!";
	open STDOUT, ">>$logfile" or die "Can't write to /dev/null: $!";
	open STDERR, ">>$logfile" or die "Can't write to /dev/null: $!";
}

sub unlock() {
	flock( LOCK, LOCK_UN ) or die "Cannot unlock - $!\n";
	close(LOCK);
	unlink($lock);
	closedir(DIR);
	close LOG;
	exit;
}

sub readlock() {
	print "Looking for LOCK at $lock\n";
	my $existing_pid;
	if ( -s $lock ) {
		open( LOCK, $lock );
		while ( my $line = <LOCK> ) {
			if ( $line =~ /(\d+)/ ) {
				my $pid_report = `ps -A |grep $1`;
				if ( $pid_report =~ /^\s*(\d+)/ ) {
					$existing_pid = $1;
				}
			}
		}
		close LOCK;
	}
	return $existing_pid;
}

sub statusme() {
	my $existing_pid = &readlock();
	if ($existing_pid) {
		print "Running as $existing_pid\n";
	} else {
		print "Found no existing program\n";
	}
}

sub killme() {
	my $existing_pid = &readlock();
	if ($existing_pid) {
		system("kill $existing_pid");
		print "Attempted to terminate existing program\n";
		return 1;
	} else {
		print "Found no existing program\n";
		return 0;
	}
}

sub process_blast($) {
	my $infile = shift;
	if ( $infile && -s $infile && -r $infile && $infile =~ /^(\w+)\.blastall.output$/) {
		my $uid = $1;
		if ($condor) {
			my $check = `grep 'Normal termination' $uid.condor.log`;
			chomp($check);
			if ( !$check ) { next OUTER; }
		} else {
			my $check = `lsof -u $user -w|grep $infile`;
			chomp($check);
			if ($check) { next OUTER; }
		}
		unless ( -s $uid . '.blastall.output.graph.png' ) {

			# make a graph
			# make csv
			# make a html
			# make tab -> bioperl broken.
			# make a gff3
			my $query_counter = 0;

			#my $writer_tab  = Bio::SearchIO::Writer::HitTableWriter->new();
			my $writer_html;
			if ($db_link) {
				$writer_html =
				  Bio::SearchIO::Writer::HTMLResultWriter->new( -nucleotide_url => $db_link,
																-protein_url    => $db_link, );
				$writer_html->id_parser( \&blast_id_parser );
			} else {
				$writer_html = Bio::SearchIO::Writer::HTMLResultWriter->new();
			}
			$writer_html->title( \&blast_title );
			$writer_html->introduction( \&blast_intro );
			$writer_html->start_report( \&blast_head );
			my $writer_gff = Bio::SearchIO::Writer::GbrowseGFF->new();
			my $writer_txt = Bio::SearchIO::Writer::TextResultWriter->new();
			my $in         = new Bio::SearchIO( -format => 'blastxml', -file => $infile );
			my $max;

			while ( my $result = $in->next_result ) {
				$max = $result->query_length()
				  if ( !$max || $max < $result->query_length() );
			}

			#rewind
			$in = new Bio::SearchIO( -format => 'blastxml', -file => $infile );
			my $graph =
			  Bio::Graphics::Panel->new(
										 -length    => $max,
										 -width     => 800,
										 -pad_left  => 10,
										 -pad_right => 10,
			  );
			while ( my $result = $in->next_result ) {
				my $query    = $result->query_name();
				my $hitcount = $result->num_hits;

				#my @hit_array;
				$query_counter++;
				if ( $hitcount && $hitcount > 0 ) {

				 #we want to know how many hits there are, only if there is a result do we bother...
				 #my $out_tab =
				 #  new Bio::SearchIO( -writer => $writer_tab,
				 #                     -file   => ">>$infile.Query$query_counter.tab" );
				 #my $out_tab_whole =
				 #  new Bio::SearchIO( -writer => $writer_tab,
				 #                     -file   => ">>$infile.tab" );
					my $out_html =
					  new Bio::SearchIO( -writer => $writer_html,
										 -file   => ">$infile.Query$query_counter.html" );
					my $out_html_whole =
					  new Bio::SearchIO( -writer => $writer_html,
										 -file   => ">>$infile.html" );
					my $out_txt =
					  new Bio::SearchIO( -writer => $writer_txt,
										 -file   => ">$infile.Query$query_counter.txt" );
					my $out_txt_whole =
					  new Bio::SearchIO( -writer => $writer_txt,
										 -file   => ">>$infile.txt" );

					#graph starts in the result section and is built through the hits...
					my $full_length =
					  Bio::SeqFeature::Generic->new(
													 -start        => 1,
													 -end          => $result->query_length,
													 -display_name => $result->query_name,
					  );
					$graph->add_track(
									   $full_length,
									   -description => 'Top 20 hits',
									   -glyph   => 'arrow',
									   -tick    => 2,
									   -fgcolor => 'black',
									   -double  => 1,
									   -label   => 1,
					);
					my $track = $graph->add_track(
						-glyph       => 'graded_segments',
						-label       => 1,
						-feature_limit => 20,
						-connector   => 'dashed',
						-bgcolor     => 'blue',
						-font2color  => 'red',
						-sort_order  => 'low_score',
						-description => sub {
							my $feature = shift;
							return unless $feature->has_tag('description');
							my ($description) = $feature->each_tag_value('description');
							my $score = $feature->score;
							"E=$score; $description";
						},
					);
					my $i = int(0);
					while ( my $hit = $result->next_hit ) {
						my $hitname = $hit->name;
						$hitname =~ s/lcl\|//;
						my $description = $hit->description;
						$description =~ s/No definition line found//;
						if ( $hit && $i < 20 ) {
							$i++;
							my $feature =
							  Bio::SeqFeature::Generic->new(
															-score        => $hit->significance,
															-display_name => $hitname,
															-tag => { description => $description },
							  );
							while ( my $hsp = $hit->next_hsp ) {
								$feature->add_sub_SeqFeature( $hsp, 'EXPAND' );
							}
							$track->add_feature($feature);
						}
					}

					# $out_tab->write_result($result);
					# $out_tab_whole->write_result($result);
					my $html_success = $out_html->write_result($result);
					$out_html_whole->write_result($result);
					$out_txt_whole->write_result($result);
					$out_txt->write_result($result);
				}
			}
			open( GRAPH, ">$infile.graph.png" );
			print GRAPH $graph->png if $graph;
			close GRAPH;
		}
	}
}

sub process_annot8r($) {
	my $software = 'annot8r';
	my $infile   = shift;
	if ( $infile =~ /^(\w+)\.$software.output$/ && -s $infile ) {
		my $uid = $1;
		if ($condor) {
			my $check = `grep 'Normal termination' $uid.condor.log`;
			chomp($check);
			if ( !$check ) { next OUTER; }
		} else {
			my $check = `lsof -u $user -w|grep $infile `;
			chomp($check);
			if ($check) { next OUTER; }
		}
		my $graph_file = $infile . ".graph.png";
		unless ( -s $graph_file ) {
			my %query_hash;
			my $in = Bio::Tools::GFF->new( -gff_version => 3, -file => $infile )
			  || die("Cannot open $infile\n");
			my $max = 100;
			while ( my $feat = $in->next_feature ) {
				unless ( $feat->source_tag eq 'ANNOT8R_physprop' ) { next; }
				$max = $feat->end() if ( $feat->end() && $feat->end() > $max );
			}
			undef($in);
			$in = Bio::Tools::GFF->new( -gff_version => 3, -file => $infile );
			my $graph =
			  Bio::Graphics::Panel->new(
										 -length    => $max,
										 -width     => 800,
										 -pad_left  => 10,
										 -pad_right => 10,
										 -key_style => 'left',
			  );
			while ( my $feat = $in->next_feature ) {
				if ( $feat->source_tag eq 'ANNOT8R_physprop' ) {
					$graph->add_track(
									   $feat,
									   -description => 'Top 20 hits',
									   -glyph   => 'arrow',
									   -tick    => 2,
									   -fgcolor => 'black',
									   -double  => 1,
									   -label   => 1,
					);
				} else {

#ID=Unknown_Query:KEGG:7;Alias=K05277;Dbxref=KEGG_ORTHOLOGY:K05277,KEGG_PATHWAY:00941,UNIPROT:Q96323;Name=Unknown_Query:KEGG:7;Note=Flavonoid biosynthesis
					my $toshow = "";
					my $source = $feat->source_tag();
					my @alias  = $feat->get_tag_values('Alias');
					my $id     = join( ' ', @alias );
					$id = 'GO terms: ' . $id   if $source eq 'ANNOT8R_GO';
					$id = 'EC terms: ' . $id   if $source eq 'ANNOT8R_EC';
					$id = 'KEGG terms: ' . $id if $source eq 'ANNOT8R_KEGG';
					$feat->seq_id($id);
					my @dbxref = $feat->get_tag_values('Dbxref');
					my @note   = $feat->get_tag_values('Note');
					$feat->add_tag_value( 'description',
										  join( ',', @dbxref ) . ' ' . join( '; ', @note ) );
					my $track = $graph->add_track(
						-glyph         => 'flag',
						-label         => 1,
						-feature_limit => 20,
						-text          => ' ',
						-category      => sub { my $f = shift; $f->source_tag() },
						-bgcolor       => sub {
							my $feature = shift;
							my $s       = $feature->source_tag();
							return 'blue'   if ( $s && $s eq 'ANNOT8R_GO' );
							return 'yellow' if ( $s && $s eq 'ANNOT8R_EC' );
							return 'green'  if ( $s && $s eq 'ANNOT8R_KEGG' );
							return 'black';
						},
						-font2color  => 'red',
						-sort_order  => 'high_score',
						-description => sub {
							my $f     = shift;
							my $score = $f->score();
							return $score unless $f->has_tag('description');
							$score .= '; '   if $score;
							$score .= 'NA; ' if !$score;
							my @desc = $f->get_tag_values('description');
							return $score . join( '; ', @desc );
						}
					);
					$track->add_feature($feat) if $feat;
				}
			}
			open( GRAPH, ">$graph_file" );
			print GRAPH $graph->png if $graph;
			close GRAPH;
			$in->close();
		}
	}
}

sub process_iprscan($) {
	my $software = 'iprscan';
	my $infile   = shift;
	if ( $infile =~ /^(\w+)\.$software.output$/ && -s $infile ) {
		my $uid = $1;
		if ($condor) {
			my $check = `grep 'Normal termination' $uid.condor.log`;
			chomp($check);
			if ( !$check ) { next OUTER; }
		} else {
			my $check = `lsof -u $user -w|grep $infile `;
			chomp($check);
			if ($check) { next OUTER; }
		}
		my $graph_file = $infile . ".graph.png";
		if ( !-s "$uid.$software.error.temp" ) {
			sleep(5);
		}
		my $job_id = `grep ^SUBMITTED $uid.$software.error.temp`;
		$job_id =~ /^SUBMITTED\s+(\S+)/;
		$job_id = $1 if $1;
		unless ( -s $graph_file ) {

			# job-wide conversions
			system( $var_hash->{'ipr_convert'}->{'exec'}
					. " -input $infile -format gff3 -output $infile.gff" );
			system( $var_hash->{'ipr_convert'}->{'exec'}
					. " -input $infile -format xml -output $infile.xml" );
			system( $var_hash->{'ipr_convert'}->{'exec'}
					. " -input $infile -format txt -output $infile.txt" );
			system( $var_hash->{'ipr_convert'}->{'exec'}
					. " -input $infile -format html -jobid $job_id > $infile.html.orig" );
			open( HTML,    "$infile.html.orig" );
			open( HTMLOUT, ">$infile.html" );
			while ( my $ln = <HTML> ) {

				if ( $ln =~ /<form method/ ) {
					while ( $ln !~ /<\/form>/ ) {
						$ln = <HTML>;
					}
					$ln = <HTML>;
				} elsif ( $ln =~ /<a href="\// ) {
					$ln =~ s/<a href="\/[^>]+>(.+)<\/a>/$1/;
					print HTMLOUT $ln;
					next;
				}
				print HTMLOUT $ln;
			}
			close(HTML);
			close(HTMLOUT);
			unlink("$infile.html.orig");
			my %query_hash;

			# query specific conversions
			open( IN, $infile );
			while ( my $ln = <IN> ) {
				my @data = split( "\t", $ln );
				push( @{ $query_hash{ $data[0] } }, $ln ) if $data[1];
			}
			close(IN);
			my $query_counter = 0;
			foreach my $id ( keys %query_hash ) {
				$query_counter++;
				open( OUT, ">$infile.Query$query_counter.raw" );
				foreach my $line ( @{ $query_hash{$id} } ) {
					print OUT $line;
				}
				close(OUT);
			}
			for ( my $i = 1 ; $i <= $query_counter ; $i++ ) {
				my $base = "$infile.Query$i";
				if ( -s $base . '.raw' ) {
					system( $var_hash->{'ipr_convert'}->{'exec'}
							. " -input $base.raw -format gff3 -output $base.gff" );
					system( $var_hash->{'ipr_convert'}->{'exec'}
							. " -input $base.raw -format xml -output $base.xml" );
					system( $var_hash->{'ipr_convert'}->{'exec'}
							. " -input $base.raw -format txt -output $base.txt" );
					system( $var_hash->{'ipr_convert'}->{'exec'}
							. " -input $base.raw -format html -jobid $job_id > $base.html.orig" );
					open( HTML,    "$base.html.orig" );
					open( HTMLOUT, ">$base.html" );
					while ( my $ln = <HTML> ) {
						if ( $ln =~ /<form method/ ) {
							while ( $ln !~ /<\/form>/ ) {
								$ln = <HTML>;
							}
							$ln = <HTML>;
						} elsif ( $ln =~ /<a href="\// ) {
							$ln =~ s/<a href="\/[^>]+>(.+)<\/a>/$1/;
							print HTMLOUT $ln;
							next;
						}
						print HTMLOUT $ln;
					}
					close(HTML);
					close(HTMLOUT);
					unlink("$base.html.orig");
				}
			}
			my $in = Bio::Tools::GFF->new( -gff_version => 3, -file => $infile . '.gff' );
			my $max = 100;
			while ( my $feat = $in->next_feature ) {
				unless ( $feat->primary_tag eq 'region' ) { next; }
				$max = $feat->end() if ( $feat->end() && $feat->end() > $max );
			}
			$in = '';
			$in = Bio::Tools::GFF->new( -gff_version => 3, -file => $infile . '.gff' );
			my $graph =
			  Bio::Graphics::Panel->new(
										 -length    => $max,
										 -width     => 800,
										 -pad_left  => 10,
										 -pad_right => 10,
			  );
			while ( my $feat = $in->next_feature ) {
				if ( $feat->primary_tag eq 'region' ) {
					$graph->add_track(
									   $feat,
									   -description => 'Top 20 hits',
									   -glyph   => 'arrow',
									   -tick    => 2,
									   -fgcolor => 'black',
									   -double  => 1,
									   -label   => 1,
					);
					my $track = $graph->add_track(
						-glyph         => 'graded_segments',
						-label         => 1,
						-feature_limit => 20,
						-connector     => 'dashed',
						-bgcolor       => sub {
							my $feature = shift;
							my $s       = $feature->score();
							return 'blue' if ( $s && $s =~ /^\d/ );
							return 'green';
						},
						-font2color  => 'red',
						-sort_order  => 'high_score',
						-description => sub {
							my $feature = shift;
							my $score   = $feature->score();
							return $score unless $feature->has_tag('description');
							$score .= '; '   if $score;
							$score .= 'NA; ' if !$score;
							my @desc = $feature->get_tag_values('description');
							return $score . join( '; ', @desc );
						}
					);
					&process_iprscan_hit($in,$track);
				}
			}
			open( GRAPH, ">$graph_file" );
			print GRAPH $graph->png if $graph;
			close GRAPH;
			$in->close();
		}
	}
}

sub process_iprscan_hit($) {
	my $feat_feed = shift                      || return;
	my $track = shift;
	while (my $hit       = $feat_feed->next_feature()){
	
	my $tag       = $hit->primary_tag();
	if ( !$tag || $tag ne 'match_set' ) { next; }

	# each iprscan hit has one match_set and one match_part after the main region file.
	$hit->seq_id( $hit->get_tag_values('Name') );
	my @array;
	my @notes = $hit->get_tag_values('Note') if $hit->has_tag('Note');
	foreach (@notes) {
		$_ =~ /(IPR\d+)/;
		$_ = $1 if $1;
	}
	push( @array, @notes ) if @notes;
	push( @array, $hit->get_tag_values('Ontology_term') )
	  if $hit->has_tag('Ontology_term');
	$hit->add_tag_value( 'description', join( '; ', @array ) ) if @array;
	my $hsp = $feat_feed->next_feature();
	$hit->add_sub_SeqFeature( $hsp, 'EXPAND' );
	my $score = $hsp->score();
	$hit->score($score) if ($score && $score ne 'NA');
        $track->add_feature($hit) if $hit;
	}
}

sub process_ssaha2($) {
	my $software = 'ssaha2';
	my $infile   = shift;
	if ( $infile =~ /^(\w+)\.$software.output$/ && -s $infile ) {
		my $uid = $1;
		if ($condor) {
			my $check = `grep 'Normal termination' $uid.condor.log`;
			chomp($check);
			if ( !$check ) { next OUTER; }
		} else {
			my $check = `lsof -u $user -w|grep $infile `;
			chomp($check);
			if ($check) { next OUTER; }
		}
		my $graph_file = $infile . ".graph.png";
		unless ( -s $graph_file ) {

			# SSAHA FORMAT
			#$VAR0 = '1279';    score
			#$VAR1 = 'query';   name of query
			#$VAR2 = '2L';      name of hit
			#$VAR3 = '1';       start of query
			#$VAR4 = '1279';    end of query
			#$VAR5 = '64311';   start of hit
			#$VAR6 = '65589';   end of hit
			#$VAR7 = 'F';       strand
			#$VAR8 = '1279';    bases matching
			#$VAR9 = '100.00';  % identity
			#$VAR10 = '1279';   length of query
			my %hash;
			my $max = 100;
			open( IN, $infile ) || warn "Cannot open $infile\n";
			while ( my $input = <IN> ) {
				unless ( $input =~ /^ALIGNMENT:/ ) { next; }
				$input =~ s/^ALIGNMENT::\d+\s+//;
				chomp($input);
				my @data = split( " ", $input );
				if ( $data[10] ) {
					$max = $data[10] if $data[10] > $max;
					push( @{ $hash{ $data[1] }{ $data[2] } }, \@data );
				}
			}
			close(IN);
			if ( !%hash ) {
				open( GRAPH, ">graph.png" );
				print GRAPH "No hits\n";
				close GRAPH;
				exit;
			}
			my $writer_html = Bio::SearchIO::Writer::HTMLResultWriter->new();
			my $writer_txt  = Bio::SearchIO::Writer::TextResultWriter->new();
			my $gff_writer =
			  Bio::Tools::GFF->new( -file        => ">>$infile.gff",
									-gff_version => 3 );
			my $out_txt_whole = ">$infile.txt";
			open( OUT_WHOLE, $out_txt_whole );
			print OUT_WHOLE
"SSAHA2:: SCORE\tQuery Name\tHit Name\tQuery Start\tQuery End\tHit Start\tHit End\tStrand\tAlignment length\tIdentity %\tQuery length\n";
			my $query_counter = int(0);
			my $graph =
			  Bio::Graphics::Panel->new(
										 -length    => $max,
										 -width     => 800,
										 -pad_left  => 10,
										 -pad_right => 10,
			  );

			foreach my $query ( keys %hash ) {
				$query_counter++;
				my $out_txt = ">$infile.Query$query_counter.txt";
				open( OUT_QUERY, $out_txt );
				print OUT_QUERY
"SSAHA2:: SCORE\tQuery Name\tHit Name\tQuery Start\tQuery End\tHit Start\tHit End\tStrand\tAlignment length\tIdentity %\tQuery length\n";
				foreach my $hit ( keys %{ $hash{$query} } ) {
					if ( $query_counter > 10 ) { last; }
					my @input = @{ $hash{$query}{$hit} };
					my $full_length =
					  Bio::SeqFeature::Generic->new(
													 -start        => 1,
													 -end          => $input[0][10],
													 -display_name => $query,
					  );
					$graph->add_track(
									   $full_length,
									   -description => 'Top 20 hits',
									   -glyph   => 'arrow',
									   -tick    => 2,
									   -fgcolor => 'black',
									   -double  => 1,
									   -label   => 1,
					);
					last;
				}
				foreach my $hit ( keys %{ $hash{$query} } ) {
					my @data   = @{ $hash{$query}{$hit} };
					my $name   = $data[0][2];
					my $strand = $data[0][7] eq 'F' ? "+1" : -1;
					my $feature =
					  Bio::SeqFeature::Generic->new(
													 -seq_id       => $query,
													 -source_tag   => 'SSAHA2',
													 -primary      => 'DNA',
													 -display_name => $name,
													 -start        => 1,
													 -end          => $data[0][10],
					  );
					$feature->add_tag_value( 'description', "Top HSP score:" . $data[0][0] );
					$feature->add_tag_value( 'description', $data[0][9] . " % identical" );
					$gff_writer->write_feature($feature);
					my $hsp_number = 0;

					foreach my $hsp (@data) {
						$hsp_number++;
						my $print = "SSAHA2:: ";
						foreach (@$hsp) { $print .= $_ . "\t"; }
						$print .= "\n";
						my $start  = @$hsp[3];
						my $end    = @$hsp[4];
						my $strand = @$hsp[7] eq 'F' ? "+1" : -1;
						if ( !$start || !$end ) { next; }
						my $subfeature =
						  Bio::SeqFeature::Generic->new(
														-display_name => $name . " hit $hsp_number",
														-seq_id       => $query,
														-source_tag   => 'SSAHA2',
														-primary      => "nucleotide_match",
														-start        => $start,
														-end          => $end,
														-strand       => $strand,
														-score        => @$hsp[0],
						  );
						$subfeature->add_tag_value( 'Target', @$hsp[2] );
						$subfeature->add_tag_value( 'Target', @$hsp[5] );
						$subfeature->add_tag_value( 'Target', @$hsp[6] );
						$subfeature->add_tag_value( 'Note',   'Ident. ' . @$hsp[9] . ' % ' );
						$feature->add_sub_SeqFeature( $subfeature, 'EXPAND' );
						$gff_writer->write_feature($subfeature);
						print OUT_WHOLE $print;
						print OUT_QUERY $print;
					}
					my $track = $graph->add_track(
						-glyph         => 'graded_segments',
						-label         => 1,
						-connector     => 'dashed',
						-bgcolor       => 'blue',
						-font2color    => 'red',
						-sort_order    => 'high_score',
						-feature_limit => 20,
						-description   => sub {
							my $feature = shift;
							return unless $feature->has_tag('description');
							return join( '; ', $feature->get_tag_values('description') );
						}
					);
					$track->add_feature($feature) if ( $query_counter <= 10 );
				}
			}
			close OUT_WHOLE;
			close OUT_QUERY;
			open( GRAPH, ">$graph_file" );
			print GRAPH $graph->png;
			close GRAPH;
		}
	}
}

sub check_submissions() {
	rewinddir(DIR);
  OUTER: while ( my $file = readdir(DIR) ) {
		if ( $file =~ /^(\w+)\.(\w+)\.wait$/ && -s $file ) {
			my $uid     = $1;
			my $softuid = $uid . '.' . $2;

			#print "Looking at file $file with $uid and $softuid\n";
			unless (
				   -s $softuid . '.done'
				|| -s $softuid
				. '.error'

				#|| -s $softuid . "*.output"
				|| -s $softuid . '.condor.submitted'
			  )
			{

				#print LOG 'Checking '.$file;
				my $software;
				my @par_array;
				my $env;
				open( FILE, $file ) || print "Cant open $file $!";
				while ( my $line = <FILE> ) {
					chomp($line);
					my @data = split( ':', $line );
					unless ( $data[1] ) {
						print "Not valid data structure\n";
						next;
					}

					# security check
					if ( $data[0] eq 'uid' && $data[1] ne $uid ) {
						open( ER, ">$softuid.error" );
						print ER "Failed: $uid ne $data[1]\n";
						close ER;
						next OUTER;
					} elsif ( $data[0] eq 'par' ) {
						push( @par_array, $data[1] );
					} elsif ( $data[0] eq 'software' ) {
						if ( !exists $software_data{ $data[1] } ) {
							print LOG "Invalid software $data[1]. Not in $software_data_str\n";
							open( ER, ">$softuid.error" );
							print ER "Failed: invalid software $data[1]\n";
							close ER;
							next OUTER;
						}
						$software = $software_data{ $data[1] }{'exec'};
					}
				}
				close(FILE);
				unless ( $software && @par_array ) {
					print LOG "No valid software or parameters found\n";
					open( ER, ">$softuid.error" );
					print ER "Failed: no valid software or parameters\n";
					close ER;
					next OUTER;
				}
				foreach my $par (@par_array) {
					my $outfile ='';
					$par =~ s/$escape_pattern//g;
					if ($par =~ /(\s\-o\w*\s+)(\S+)/){
					   $outfile = $2;
					}
					if ( !$outfile ) {
						$par =~ /\s\>\s+(\S+)/;
						$outfile = $1;
					}
					my $log = "Starting $software with $par to produce output $outfile ";
					if ( !$outfile ) {
						open( ER, ">$softuid.error" );
						print $log . "\nStopped because cannot determine outfile from:\n $par.\n";
						print ER $log
						  . "\nStopped because cannot determine outfile from:\n $par.\n";
						close ER;
						next;
					}
					my $delete_file;
					if ( $software =~ /ssaha/ ) {
						$par =~ /-save\s+(\S+)/;
						print "$1 is from $par\n";
						my $subject_db = $1;
						if ($subject_db) {
							if ( !-s $subject_db . ".head" && -s $subject_db . ".head.gz" ) {

								#testing write
								my $dirname = `dirname $subject_db`;
								chomp($dirname);
								unless ( -w $dirname ) {
									open( ER, ">$softuid.error" );
									print "Cannot uncompress into $dirname: Skipping\n";
									print ER "Cannot uncompress into $dirname: Skipping\n";
									close ER;
									next;
								}
								system("gunzip -c $subject_db.head.gz >$subject_db.head");
								$delete_file = "$subject_db.head";
							}
						}
					}
					my $failed;
					if ($condor) {
						my $me = `id`;
						if ( !-s "$softuid.condor.submitted" ) {
							$log .= " using condor\n";
							print LOG $log;
							&prepare_condor( $softuid, $software, $par );
							if ( -s "$softuid.condor" ) {
								system("$condor $softuid.condor 1>> $logfile 2>>$logfile");
								system("date >$softuid.condor.submitted");
							} else {
								open( ERROR, ">$softuid.error" );
								print ERROR "Failed to find condor parameter file\n";
								close ERROR;
								$failed = 1;
							}
						}
					} else {
						$log .= " without condor\n";
						print LOG $log;
						open( ERROR, ">$softuid.error" );
						if ( $software =~ /blast/ ) {
							system("$software $par -a $cpus >>$logfile 2>>$softuid.error.temp");
						} elsif ( $software =~ /ssaha/ ) {

							#ssaha output to STDOUT so no logfile
							system("$software $par 2>>$softuid.error.temp");
						} else {
							system("$software $par >>$logfile 2>>$softuid.error.temp");
						}
						print LOG " Task completed as $outfile.\n";
						if ( $? == -1 ) {
							print ERROR "failed to execute: $!\n";
						} elsif ( $? & 127 ) {
							printf ERROR
							  "child died with signal %d, %s coredump\n",
							  ( $? & 127 ), ( $? & 128 ) ? 'with' : 'without';
						} elsif ( $? != 0 ) {
							printf ERROR "child exited with value %d\n", $? >> 8;
						}
					}
					if ( -s $softuid . '.error.temp' ) {
						open( TEMPERR, $softuid . '.error.temp' );
						while ( my $ln = <TEMPERR> ) {
							if (    $ln =~ /FATAL ERROR/i
								 || $ln =~ /had length 0/
								 || $ln =~ /^Cannot/ 
								 || $ln =~ /Some jobs failed/
								 )
							{
								$failed = 1;
								last;
							}
						}
						close(TEMPERR);
						system("mv -f $softuid.error.temp $softuid.error")
						  if $failed;
					} else {
						unlink( $softuid . '.error.temp' );
					}
					close(ERROR);
                    sleep(2);
					unless ( -s $softuid . '.error' ) {
						my $date = `date >> $softuid.wait`;
						rename( $softuid . '.wait', $softuid . '.done' );
						unlink("$softuid.error");
					}

					#delete $delete_file if it has not been accessed for 7 days
					if ( $delete_file && -s $delete_file && -A $delete_file > 7 ){
					   unlink($delete_file);
					}
				}    #end par
			}
		}
	}
}

sub check_results() {
	rewinddir(DIR);
  OUTER: while ( my $outfile = readdir(DIR) ) {
		if ( $outfile =~ /\.output$/ && -s $outfile ) {
			&process_blast($outfile)   if $var_hash->{'blastall'}->{'exec'};
			&process_ssaha2($outfile)  if $var_hash->{'ssaha2'}->{'exec'};
			&process_iprscan($outfile) if $var_hash->{'ipr_convert'}->{'exec'};
			&process_annot8r($outfile) if $var_hash->{'annot8r'}->{'exec'};
		}
	}
}

sub change_uid() {
	my $user = shift;
	my $me   = `whoami`;
	chomp($me);
	my $uid = `id -u $user 2>/dev/null`;
	chomp($uid);
	my $gid = `id -g $user 2>/dev/null`;
	chomp($gid);
	if ( !$uid || !$gid ) {
		die "Cannot find user $user\n";
	}
	if ( $me ne $user ) {
		## repeat group in effective gid to empty groups - Otherwise condor does not work
		$) = "$gid $gid";
		$( = $gid;
		setuid($uid);
		if ( $> ne $uid ) { die "Failed to change user to $uid, exiting.\n" }
	}
	print "Will run as $user ($uid,$gid)\n";
	sleep(1);
	return ( $uid, $gid );
}

sub get_uid($) {
	my $user = shift;
	my $uid  = `id -u $user`;
	chomp($uid);
	my $gid = `id -g $user`;
	chomp($gid);
	return ( $uid, $gid );
}

sub connecttodb() {
	if ( $dsn && -f $dsn ) {
		open( IN, $dsn );
		$dsn = "";
		while ( my $line = <IN> ) {
			unless ( $line =~ /^\s*$/ ) { chomp($line); $dsn .= $line; }
		}
		close(IN);
	}
	if ( !$dsn ) {
		die "No connection to database provided.\n";
	}
	my $dbh = DBI->connect($dsn);
	die("Error: Unable to connect to database") if ( !$dbh );
	return ($dbh);
}

sub disconnectdb ($) {
	my $dbh = shift;
	$dbh->disconnect();
}

sub get_pars() {
	if ( $dsn && -f $dsn ) {
		open( IN, $dsn );
		$dsn = "";
		while ( my $line = <IN> ) {
			unless ( $line =~ /^\s*$/ ) { chomp($line); $dsn .= $line; }
		}
		close(IN);
	}
	if ( !$dsn ) {
		die "No connection to database provided.\n";
	}
	my $dbh = &connecttodb($dsn)
	  || die("Cannot connect to Drupal database...\n");

	#get condor
	my $select_condor_sql =
	  "SELECT active from " . $dbprefix . "gmod_dbsf_software where uniquename='condor'";
	my $select_condor_prepare = $dbh->prepare($select_condor_sql);
	$select_condor_prepare->execute();
	my $condor_hash = $select_condor_prepare->fetchrow_hashref;
	if ( $condor_hash->{'active'} == 1 ) {
		print
"Utilizing Condor Distributed Computing service: Allowing all-write permissions for $dir\n";
		chmod( 0777, $dir );
		$condor = `which condor_submit`;
		chomp($condor);
		if ( !$condor ) {
			die
"I'm sorry, but your condor has been set to active but cannot find the environmental variables for Condor. Are you sure the Apache user has sourced the condor env. variables? Also please don't run this using sudo but if needed, login as root.\n";
		}
	}
	$select_condor_prepare->finish();

	# get executables
	my $select_software_sql =
	  "SELECT uniquename from " . $dbprefix . "gmod_dbsf_software where active is true";
	my $select_software_prepare = $dbh->prepare($select_software_sql);
	$select_software_prepare->execute();
	while ( my $results_hash = $select_software_prepare->fetchrow_hashref ) {
		$software_data{ $results_hash->{'uniquename'} }{'exists'} = 1;
		$software_data_str .= $results_hash->{'uniquename'} . ',';
	}
	die "No software have been activated using the administrator interface!\n" if !$software_data_str;
	chop($software_data_str);
	$select_software_prepare->finish();
	my $select_software_setting_sql =
	    "SELECT value FROM "
	  . $dbprefix
	  . "gmod_dbsf_softwareprop where software_id="
	  . "(SELECT software_id from "
	  . $dbprefix
	  . "gmod_dbsf_software where uniquename=?)"
	  . ' AND rank=0 AND type_id='
	  . "(SELECT cvterm_id from "
	  . $dbprefix
	  . "gmod_dbsf_cvterm as cvterm JOIN "
	  . $dbprefix
	  . "gmod_dbsf_cv as cv on cv.cv_id=cvterm.cv_id"
	  . " where cv.name='software_setting' AND cvterm.name=?)";
	my $select_software_setting_prepare = $dbh->prepare($select_software_setting_sql);
	foreach my $software ( keys %software_data ) {
		$select_software_setting_prepare->execute( $software, 'executable' );
		my $results_hash1 = $select_software_setting_prepare->fetchrow_hashref;
		if ( $results_hash1->{'value'} ) {
			$software_data{$software}{'exec'} = $results_hash1->{'value'};
		} else {
			die "Your database does not provide an executable for $software\n";
		}
		delete( $software_data{$software}{'exists'} );
		$select_software_setting_prepare->execute( $software, 'data' );
		my $results_hash2 = $select_software_setting_prepare->fetchrow_hashref;
		$software_data{$software}{'env'} = $results_hash2->{'value'}
		  if $results_hash2->{'value'};
		if ( $software eq 'blastall' ) {
			$ENV{'BLASTMAT'} = $software_data{$software}{'env'}
			  if $software_data{$software}{'env'};
		}
	}
	$select_software_setting_prepare->finish();
	&disconnectdb($dbh);
	if ($software_data_str) {
		print "Available software: $software_data_str\n";
	} else {
		die "Fatal: will not start: no software found.\n";
	}

	#just in case
	system("touch $dir/index.html");
	print "The following data is stored in your database:\n";
	foreach my $software ( keys %software_data ) {
		print " $software:\t" . $software_data{$software}{'exec'} . "\n";
	}
	return \%software_data;
}

sub blast_title($) {
	my $data = shift;
	my $alg  = $data->{'_algorithm_version'};
	return "<b>$alg</b><h3>BLAST results</h3>";
}

sub blast_head($) {
	return '
    <HTML>
      <HEAD> 
      <TITLE>BLAST results</TITLE>
      </HEAD>
      <!------------------------------------------------------------------->
      <!-- Generated by Bio::SearchIO::Writer::HTMLResultWriter          -->
      <!-- http://bioperl.org                                            -->
      <!------------------------------------------------------------------->
      <BODY BGCOLOR="WHITE">
    ';
}

sub blast_id_parser($) {
	my ($string) = @_;
	my ( $gi, $acc );
	if ( $string =~ s/gi\|(\d+)\|?// ) {
		$gi  = $1;
		$acc = $1;
	} elsif ( $string =~ /lcl\|([\w\.\_\-]+)/ ) {
		$acc = $1;
	} elsif ( $string =~ /gnl\|[\w\.\_\-]+\|([\w\.\_\-]+)/ ) {
		$acc = $1;
	} elsif ( $string =~ /(\w+)\|([\w\.\_]+)\|[A-Z\d\_]+?/ ) {
		$acc = defined $2 ? $2 : $1;
	} else {
		$acc = $string;
		$acc =~ s/^\s+(\S+)/$1/;
		$acc =~ s/\s+$//;
	}
	return ( $gi, $acc );
}

sub blast_intro($) {
	my $data = shift;
	my $db   = basename( $data->database_name );
	$db =~ s/\w+\.subject/user_uploaded_database/g;
	return '<p>' . sprintf(
		qq{
    <b>Query=</b> %s %s<br><dd>(%s letters)</dd>
    <p>
    <b>Database:</b> %s<br><dd>%s sequences; %s total letters<p></dd>
    <p>
    }, $data->query_name,
		$data->query_description,
		&_numwithcommas( $data->query_length ),
		$db,
		&_numwithcommas( $data->database_entries ),
		&_numwithcommas( $data->database_letters ),
	) . '</p>';
}

sub _numwithcommas($) {

	# from Perl Cookbook 2.17
	my $num = reverse( $_[0] );
	$num =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $num;
}

sub get_chado_seq_data($$) {
	my $dbh = &connecttodb($dsn) || die("Cannot connect to Drupal database...\n");
	my $feature_id = shift || return;
	my $feature_select_sql =
	  "SELECT uniquename,residues from " . $dbprefix . "gmod_dbsf_feature where feature_id=?";
	my $select_prepare = $dbh->prepare($feature_select_sql);
	$select_prepare->execute($feature_id);
	my $results_ref = $select_prepare->fetchrow_hashref;
	my $seq         = $results_ref->{'residues'};
	my $id          = $results_ref->{'uniquename'};
	$select_prepare->finish();
	return ( $id, $seq );
}

sub prepare_condor($$$) {
	my $uid  = shift || warn("Cannot find UID\n");
	my $exec = shift || warn("Cannot find executable\n");
	my $data = shift || warn("No data");
	if ( !$uid || !$exec || !$data ) { return; }
	my @dbs = glob("$uid.subject*");
	my $tr_str = join( ",", @dbs ) if @dbs;
	$tr_str .= ",$uid.query";
	open( CONDOR, ">$uid.condor" ) || warn "Cannot not create $uid.condor $!\n";
	print CONDOR "Universe = vanilla\n";
	print CONDOR "Executable = $exec\n";
	print CONDOR "Initial_dir = $dir\n";
	print CONDOR "Notification = NEVER\n";
	print CONDOR "should_transfer_files = NO\n";

	#    print CONDOR "should_transfer_files = YES\n";
	#    print CONDOR "transfer_input_files = $tr_str\n";
	#    print CONDOR "transfer_output_files = $uid.output\n";
	#    print CONDOR "when_to_transfer_output = ON_EXIT\n";
	print CONDOR "transfer_executable = FALSE\n";
	print CONDOR "log = $uid.condor.log\n";

	#condor is absurd with error file permissions. print CONDOR "error = $uid.condor.err\n";
	print CONDOR "getenv = True\n";
	print CONDOR "Arguments = \"$data\"\n";
	print CONDOR " Queue\n";
	close CONDOR;
	## dealing with condor's peculiarities...
	#system ("touch $dir/$uid.condor.err");
	#chmod (0666,"$dir/$uid.condor.err");
	return 1;
}

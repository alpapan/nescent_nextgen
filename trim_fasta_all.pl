#!/usr/bin/perl -w
package trim_fasta_all;
use strict;
use Data::Dumper;
our $VERSION = '1.0';

#04MAR11: Added GC/AT ratio check as ratio cutoff

=head1 NAME

trim_fasta_all.pl - removes sequences from a FASTA file 

=head1 VERSION

 Version 0.2

=head1 SYNOPSIS

trim_fasta_all.pl [options] <infiles>

removes sequences from a FASTA file. See perldoc for more info.

	'i|fa|fasta=s'    => FASTA file to trim. You can also give multiples as arguments without any -i/-fa option.
	'outfile:s'	=> Optionally, the name of the trimmed outfile
	'blastfile:s'	=> BLASTFILE to retrieve sequences from
	'blastquery'		=> grab BLAST queries 
	'blasthit'       =>  grab BLAST hits
	'evalue=s'	=> Evalue cut-off for blastfile (currently broken)
	'c|character=s' => Characters to look for. If present, remove sequence.
	'le|length=i'   => Number of minimum characters
	'p|proportion'  => Discard sequences for which a mononucleotide frequency exceeds this proportion 
	'ratio'		=> Discard sequences for which the GC or AT frequency exceeds this ratio
	'x'             => Do not include the Xx characters when calculating size of sequence
	'npl'           => Do not include these characters when calculating size: NPLnpl
	'lc|lowercase'  => Do not include lowercase characters when calculating size of sequence (e.g. to not include low quality bases)
	'id|idfile=s'   => A second FASTA file containing IDs to remove from FASTA file. Alternatively a text file with one ID per line
	'descr'		=> For above: search description line instead of primary id.
	'ci'		=> Case insensitivity for above two options
	'invert'	=> Invert match (invert output filenames)
	'log'           => Keep a log file
	'df'            => Do not write discarded sequences (less IO)
    'solq'          => Input is FASTQ (Solexa 1.3-1.4)
    'sanq'          => Input is FASTQ (Sanger)
    'single'	    => Entire output sequence/quality is in a single line (good for parsing)

=head1 DESCRIPTION

Processes file (-fa) when certain character(s) are present (-c); or a list of IDs is provided (-id); or a certain length-cut off is not satisfied (-le); or a proportion of nucleotide frequence can be specified (-p) instead. The -log option produces a log file reporting what happened to each sequence 
The option to not include Xs and/or NPLs and/or lower-case characters in the cut-off calculation is forced with -x and/or -npl and/or -lc respectively.
Uses BioPerl. A disk-friendly function (-df) prevents the FASTA file of discarded sequences of being written.

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

None known so far.

=cut
use Bio::SeqIO;
use Bio::SearchIO;
use Time::Progress;
use Getopt::Long;
use Pod::Usage;
$| = 1;
my (
	 $character,    @infiles,     $length_cutoff, $xmask,
	 $nplmask,      $ci,          $blastfile,     $evalue_cutoff,
	 $lcmask,       $prop_cutoff, @idfiles,       $log,
	 $logfile,      $invert,      $sangerfastq,   $blast_hit,$blast_query,
	 $user_outfile, $df,          %ids,           $help,
	 $convert2uc,   $descr_flag,  $solexafastq,   $search_accession,
	 $seq_search,   $single_line, $ratio_cutoff
);
GetOptions(
	'i|fa|fasta=s{,}' => \@infiles,
	'blastfile=s'     => \$blastfile,
	#'evalue=s'	=> \$evalue_cutoff,
	'c|character=s'   => \$character,
	'le|length=i'     => \$length_cutoff,
	'p|proportion=f'  => \$prop_cutoff,
	'ratio=f'	  => \$ratio_cutoff,
	'x'               => \$xmask,
	'uc|uppercase'    => \$convert2uc,
	'npl'             => \$nplmask,
	'lc|lowercase'    => \$lcmask,
	'ids|idfile=s{,}' => \@idfiles,
	'descr'           => \$descr_flag,
	'invert'          => \$invert,
	'ci'              => \$ci,
	'log'             => \$log,
	'df'              => \$df,
	'h|help'          => \$help,
	'solq'            => \$solexafastq,
	'sanq'            => \$sangerfastq,
	'seq'             => \$seq_search,
	'outfile:s'       => \$user_outfile,
	'single'          => \$single_line,
	'blastquery'     => \$blast_query,
	'blasthit'       => \$blast_hit,

	  #'accessions'=> \$search_accession,
);
if ($help) { pod2usage; }
@infiles = @ARGV if !@infiles;
unless (@infiles) {
	print "Failed to provide or find input file\n";
	pod2usage;
}
unless (    $character
		 || $length_cutoff
		 || $prop_cutoff
		 || @idfiles
		 || $blastfile )
{
	die("Nothing to do!\n");
}
unless ($evalue_cutoff) { $evalue_cutoff = 1; }
my $timer   = new Time::Progress;
my $counter = int(0);
foreach my $idfile (@idfiles) {
	if ( $idfile && -s $idfile ) {
		my $pattern;
		if   ($descr_flag) { $pattern = '^\s*\S+\s+(.+)$'; }
		else               { $pattern = '^\s*(\S+)\s*'; }
		my @test_lines = `head $idfile`;
		foreach my $test (@test_lines) {
			if ( $test =~ /^>/ ) { $pattern = "Bio::SeqIO"; }
		}
		my $number = `wc -l < $idfile`;
		chomp($number);
		$number /= 2 if $pattern eq "Bio::SeqIO";
		$timer->attr( min => 0, max => $number );
		$timer->restart;
		print "Building hash from $idfile with $pattern ($number lines)\n";
		my $flag;

		if ( $pattern eq "Bio::SeqIO" ) {
			my $id_obj = new Bio::SeqIO( -file => $idfile, -format => "fasta" );
			while ( my $object = $id_obj->next_seq() ) {
				$counter++;
				print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /00000$/ );
				if    ($seq_search) { $ids{ $object->seq() }         = 1; }
				elsif ($descr_flag) { $ids{ $object->description() } = 1; }
				else                { $ids{ $object->id() }          = 1; }
				$flag = 1 if !$flag;
			}
		} else {
			open( IN, $idfile ) || die();
			while ( my $line = <IN> ) {
				$counter++;
				print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /00000$/ );
				if ($ci) {
					if ( $line =~ /$pattern/i ) {
						$ids{$1} = 1;
						$flag = 1 if !$flag;
					}
				} else {
					if ( $line =~ /$pattern/ ) {
						$ids{$1} = 1;
						$flag = 1 if !$flag;
					}
				}
			}
			close(IN);
		}
		if ( !$flag ) { die "Failed to get list of IDs to extract...\n"; }
		else {
			print "Hash presence of $idfile verified\n";
		}
	} elsif ($idfile) {
		warn "File $idfile is empty or does not exist!\n";
	}
}
if ( $blastfile && -s $blastfile ) {
    if ($blast_hit){
	 print "Building HASH for queries and hits from $blastfile...\n";
	 my @blast_hits = `grep '^>' $blastfile`;
      $timer->attr( min => 0, max => scalar(@blast_hits) );
      $timer->restart;
      chomp(@blast_hits);
      foreach my $blast (@blast_hits) {
      #next if $blast=~/^Sbjct|^Query|^Number|^Matrix:|^Gap penalties|^Length|^Database|^BLASTN|^Jinghui|^Database|^programs/i;
        $counter++;
        if ( $counter =~ /0000$/ ) {
            print $timer->report( "eta: %E min, %40b %p\r", $counter );
        }
        $blast =~ /^>(\S+)/;
        $ids{$1} = 1;
      }
      print "Found $counter significant results\n";
    }elsif($blast_query){
	print "Building HASH for queries from $blastfile...\n";
	my @blast_queries = `grep -B 18 '^Sequences producing' $blastfile |grep '^Query='`;
	
	$timer->attr( min => 0, max => scalar(@blast_queries) );
	$timer->restart;
	chomp(@blast_queries);
	foreach (@blast_queries) {
	  next if $_=~/^Sbjct|^Query|^Number|^Matrix:|^Gap penalties|^Length/i;
		$counter++;
		if ( $counter =~ /0000$/ ) {
			print $timer->report( "eta: %E min, %40b %p\r", $counter );
		}
		$_ =~ s/^Query=\s+//;
		$ids{$_} = 1;
	}
	print "Found $counter significant results\n";
    }else{
      die "Please provide -blasthit and/or -blastquery\n";
    }
}
foreach my $file (@infiles) {
	&process($file);
}
#####################################################
sub process ($) {
	my $fastafile = shift;
	print "Processing... $fastafile\n";
	my $fastafiletrim = "$fastafile.trim";
	$fastafiletrim = $user_outfile if $user_outfile;
	my $fastafilediscard = "$fastafile.discard";
	$fastafilediscard = $user_outfile . ".discard" if $user_outfile;
	if (!-s $fastafile){
		warn "File not found, skipping\n";
		return;
	}if (-s $fastafiletrim){
		warn "Output file $fastafiletrim already exists, skipping\n";
		return;
	}
	my ( $filein, $fileout, $fileout2, $number );
	if ($solexafastq) {
		$number = `grep -c '^\@' $fastafile`;
		$filein =
		  new Bio::SeqIO( -file => $fastafile, -format => "fastq-solexa" )
		  unless $single_line;
		$fileout =
		  new Bio::SeqIO( -file => ">$fastafiletrim",
						  -format => "fastq-solexa" )
		  unless $single_line;
		$fileout2 = new Bio::SeqIO(
									-file   => ">$fastafilediscard",
									-format => "fastq-solexa"
		) unless $single_line;
		open( IN,   $fastafile )           if $single_line;
		open( OUT1, ">$fastafiletrim" )    if $single_line;
		open( OUT2, ">$fastafilediscard" ) if $single_line;
	} elsif ($sangerfastq) {
		$number = `grep -c '^\@' $fastafile`;
		$filein = new Bio::SeqIO( -file => $fastafile, -format => "fastq" )
		  unless $single_line;
		$fileout =
		  new Bio::SeqIO( -file => ">$fastafiletrim", -format => "fastq" )
		  unless $single_line;
		$fileout2 =
		  new Bio::SeqIO( -file => ">$fastafilediscard", -format => "fastq" )
		  unless $single_line;
		open( IN,   $fastafile )           if $single_line;
		open( OUT1, ">$fastafiletrim" )    if $single_line;
		open( OUT2, ">$fastafilediscard" ) if $single_line;
	} else {
		$number = `grep -Fc '>' $fastafile`;
		$filein = new Bio::SeqIO( -file => $fastafile, -format => "fasta" )
		  unless $single_line;
		$fileout =
		  new Bio::SeqIO( -file => ">$fastafiletrim", -format => "fasta" )
		  unless $single_line;
		$fileout2 =
		  new Bio::SeqIO( -file => ">$fastafilediscard", -format => "fasta" )
		  unless $single_line;
		open( IN,   $fastafile )           if $single_line;
		open( OUT1, ">$fastafiletrim" )    if $single_line;
		open( OUT2, ">$fastafilediscard" ) if $single_line;
	}
	chomp($number);
	if ($log) {
		$logfile = $fastafile . ".trim.log";
		open( LOG, ">$logfile" );
	}
	if ($log) { print LOG "FASTA $fastafile contained $number sequences\n"; }
	my ( $empty, $discard, $trim );
	$counter = 0;
	$timer   = new Time::Progress;
	$timer->attr( min => 0, max => $number );
	$timer->restart;
	while ( my $object = $single_line ? <IN> : $filein->next_seq() ) {
	        next if !$object;
       		next if $single_line && $object=~/^\s*$/;
		$counter++;
		if ( $counter =~ /0000$/ ) {
			print $timer->report( "eta: %E min, %40b %p\r", $counter );
		}
		my ( $id, $sequence, $description, $qual );
		if ($single_line) {
			chomp($object);
			$object =~ /^(\S)(\S+)\s*(\S*)/;
			$id          = $2;
			$description = $3;
			$sequence    = <IN>;
			chomp($sequence);
			die "Sequence $counter has a header which starts with $1. This does not seem to be right...\n$object\n$sequence\n" unless ($1 eq '>' || $1 eq '@' || $1 eq '+');
			if ( $solexafastq || $sangerfastq ) {
				$qual = <IN> . <IN>;
				chomp($qual);
			}
		} else {
			$id          = $object->id();
			$sequence    = $object->seq() if ($seq_search);
			$description = $object->description() ? $object->description() : '';
		}

		# trim if given an ID file
		if ( @idfiles || $blastfile ) {
			if ( $sequence && $seq_search ) {
				if ( $ids{$sequence} ) {
						unless ( $df && !$invert ) {
							if ($single_line) {
								if ($qual) {
									print OUT2 "@" . "$id\n$sequence\n$qual\n";
								} else {
									print OUT2 ">$id\n$sequence\n";
								}
							} else {
								$fileout2->write_seq($object);
							}
						}
						$discard++;
						if ($log) {
							print LOG "Sequence $id discarded because the Sequence was found in idfiles\n";
						}
						#DO get it more than once
						#delete($ids{$sequence});
						next;
				} else {
					next;
				}
			} elsif ( exists $ids{$id} && $ids{$id}==1) {
				unless ( $df && !$invert ) {
					if ($single_line) {
						if ($qual) {
							print OUT2 "@" . "$id\n$sequence\n$qual\n";
						} else {
							print OUT2 ">$id\n$sequence\n";
						}
					} else {
						$fileout2->write_seq($object);
					}
				}
				$discard++;
				if ($log) {
					print LOG "Sequence $id discarded because the ID was found in idfiles\n";
				}

				#make sure we don't get it twice
				$ids{$id}=2;
				next;
		    } elsif ( exists $ids{$id}) {
		    	next;
				# if id exists multiple times don't write it in any file.
			} elsif ( exists $ids{ $id . ' ' . $description } && $ids{ $id . ' ' . $description }==1) {
				unless ( $df && !$invert ) {
					if ($single_line) {
						if ($qual) {
							print OUT2 "@" 
							  . $id
							  . $description
							  . "\n$sequence\n$qual\n";
						} else {
							print OUT2 ">" 
							  . $id
							  . $description
							  . "\n$sequence\n";
						}
					} else {
						$fileout2->write_seq($object);
					}
				}
				$discard++;
				if ($log) {
					print LOG "Sequence $id.$description discarded because the ID was found in idfiles\n";
				}

				#make sure we don't get it twice
				$ids{ $id . ' ' . $description } =2;
				next;
			} elsif ( exists $ids{ $id . ' ' . $description }) {
				 next;
			}
		}
		$sequence = $object->seq() if !$sequence;
		if ($sequence) {
			if ($xmask)   { $sequence =~ s/[X]//ig; }
			if ($nplmask) { $sequence =~ s/[NPL]//ig; }
			if ($lcmask)  { $sequence =~ s/[a-z]//g; }
			my $length = length($sequence);

			# trim if given a character(s)
			if ($character) {
				if ( $sequence =~ /[$character]/ ) {
					unless ( $df && !$invert ) {
						if ($single_line) {
							if ($qual) {
								print OUT2 "@" . "$id\n$sequence\n$qual\n";
							} else {
								print OUT2 ">$id\n$sequence\n";
							}
						} else {
							$fileout2->write_seq($object);
						}
					}
					$discard++;
					if ($log) {
						print LOG
"Sequence $id discarded because character $character was found\n";
					}
					next;
				}
			}

			#trim if given a length cutoff
			if ($length_cutoff) {
				if ( !$length || $length < $length_cutoff ) {
					unless ( $df && !$invert ) {
						if ($single_line) {
							if ($qual) {
								print OUT2 "@" . "$id\n$sequence\n$qual\n";
							} else {
								print OUT2 ">$id\n$sequence\n";
							}
						} else {
							$fileout2->write_seq($object);
						}
					}
					$discard++;
					if ($log) {
						print LOG
"Sequence $id discarded because length $length was smaller than cutoff $length_cutoff\n";
					}
					next;
				}
			}

			#trim if given a proportion of A/T/C/G
			if ($prop_cutoff || $ratio_cutoff) {
				my $As = ( $sequence =~ tr/A// );
				my $Ts = ( $sequence =~ tr/T// );
				my $Cs = ( $sequence =~ tr/C// );
				my $Gs = ( $sequence =~ tr/G// );
				my $propA   = ( $As / $length );
				my $propT   = ( $Ts / $length );
				my $propC   = ( $Cs / $length );
				my $propG   = ( $Gs / $length );
				my $GCratio = $propG + $propC if $ratio_cutoff;
				my $ATratio = 1 - $GCratio if $ratio_cutoff;
				if (  $prop_cutoff &&( 
					 $propA > $prop_cutoff
					 || $propT > $prop_cutoff
					 || $propG > $prop_cutoff
					 || $propC > $prop_cutoff )
				   || $ratio_cutoff && (
					 $ATratio > $ratio_cutoff
					 || $GCratio > $ratio_cutoff )
				   )
				{

					unless ( $df && !$invert ) {
						if ($single_line) {
							if ($qual) {
								print OUT2 "@" . "$id\n$sequence\n$qual\n";
							} else {
								print OUT2 ">$id\n$sequence\n";
							}
						} else {
							$fileout2->write_seq($object);
						}
					}
					$discard++;
					if ($log) {
						print LOG "Sequence $id discarded because of one nucleotide proportion (A:$propA; T:$propT; G:$propG; C:$propC higher than cutoff $prop_cutoff or GC/AT higher than $ratio_cutoff\n" if $ratio_cutoff && $prop_cutoff;
						print LOG "Sequence $id discarded because of GC/AT proportion (A:$propA; T:$propT; G:$propG; C:$propC) higher than $ratio_cutoff\n" if $ratio_cutoff;
						print LOG "Sequence $id discarded because of one nucleotide proportion (A:$propA; T:$propT; G:$propG; C:$propC higher than cutoff $prop_cutoff\n" if $prop_cutoff;
					}
					next;
				}
			}

			#next has taken care of discards.
			$trim++;
			if ($convert2uc) {
				$object->seq( uc($sequence) ) if !$single_line;
				$sequence = uc($sequence) if $single_line;
			}
			unless ( $df && $invert ) {
				if ($single_line) {
					if ($qual) {
						print OUT1 "@" . "$id\n$sequence\n$qual\n";
					} else {
						print OUT1 ">$id\n$sequence\n";
					}
				} else {
					$fileout->write_seq($object);
				}
			}
		}    #end if $sequence
		else {
			$empty++;
			if ($log) {
				print LOG "Sequence $id discard because it was empty\n";
			}
			next;
		}
	}
	if ( !$empty )   { $empty   = int(0); }
	if ( !$discard ) { $discard = int(0); }
	if ( !$trim )    { $trim    = int(0); }
	if ($invert) {
		system("mv -i $fastafiletrim tmpfile");
		system("mv $fastafilediscard $fastafiletrim");
		system("mv tmpfile $fastafilediscard");
		my $temp = $trim;
		$trim    = $discard;
		$discard = $temp;
	}
	unless ( -s "$fastafilediscard" ) { unlink "$fastafilediscard"; }
	my $elapsed = $timer->report("%L");
	print
"\nTime elapsed: $elapsed min.\nDone, $empty were empty and an additional $discard were discarded. Kept $trim as $fastafiletrim\n";
	if ($log) {
		print LOG
"\n$empty were empty and an additional $discard were discarded. Kept $trim as $fastafiletrim\n";
	}
	close(LOG);
}
print "\n";

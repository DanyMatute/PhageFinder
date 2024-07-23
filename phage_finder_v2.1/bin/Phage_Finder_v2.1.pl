#!/usr/bin/env perl
#
#Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ver 2.0 to fix failure to properly loop through every contig.  Added HMMER3 input capability. Standardized region .pep, .seq, and .con output filenames.  Changed .tab to PFPR_tab.txt to easy opening in MS Excel. Added curated lists of core phage functions to aid in identification of these genes.
# ver 1.9beta to fix # of hits per window that is allowable to be considered as a valid phage region (gdv phage region 8 had 3 hits per window)
# ver 1.8 12/30/02: added capability to search for phage HMM hits and include data in output files
# ver 1.8 12/18/02: restructured subs find_5prime_end and find_3prime_end to look for hits slightly upstream of orf (within window or next window) and use an array to move through end5s of orfs
# ver 1.7 11/14/02: fix problem where 1 window with no BLAST hits within a larger phage region prematurely terminates the region - look at next window as well (sub find_regions) 
# ver 1.6 11/14/02: check for existence of ($phage_hash{1}->{'seed_end5'}) [ie a putative hit region] so that the program will not waist time extending regions that do not exist
# ver 1.6 streamlined counting loop within MAIN to delete hash keys of $searchhash when already counted and also end looping through keys of searchhash when keys >= window size
# ver 1.6 fixed problem where only one role_id was stored per ORF (now $hithash and $rehash can hold multiple role_ids)
# ver 1.2 fills in all featnames into %hithash (in sub select_featnames_from_btab), but only hit=1 if $a[10] != -100
# ver 1.2 use SQL to only pull end5s of orfs in specific role categories (typical of phage regions) and record that role_id
# ver 1.3 will add a max peak height counter and store median of max peak in $fighash{$n}->{'median'} (will later be used for outward scans to find more potential phage genes
# ver 1.1 includes step size separate from window size (like selfsim)
#
# Boolean logic convention True = 1, False = 0

my $prog = $0;
$prog =~ s/.*\///;
my $invocation = $prog . " @ARGV";
my $version;
use strict;
BEGIN { # set version variable before use lib at compile time
    $version = "2.1";
}
use lib "$ENV{'HOME'}/phage_finder_v$version/lib/";
use Math::Round qw(:all);
use PHAGE::Phage_subs_v2;
use Getopt::Std;
getopts ('hdVb:p:t:w:s:P:U:E:H:a:A:B:i:m:r:n:F:S');

our ($opt_h,$opt_b,$opt_p,$opt_t,$opt_i,$opt_m,$opt_w,$opt_s,$opt_E,$opt_H,$opt_a,$opt_r,$opt_n,$opt_B,$opt_A,$opt_F,$opt_S,$opt_d,$opt_V);
my ($basedir,$default_hpw,$defined_asmbl_id,$window,$step,$evalue,$asmbl_id,$btabpath,$btabfile,$infofile,$hmmfile,$tRNA_data,$tmRNA_data,$asmbly_file_path,$asmbly_file,$search_method,$strict,$project_prefix,$DEBUG);

if ($opt_h) { &option_help; }
if ($opt_V) {die "$version\n";}
if (length($opt_b) >0) {$basedir = $opt_b; print "$basedir\n";} else { $basedir = $ENV{'PWD'};} # if no value for option b (base or working directory), quit with help menu
if (length($opt_p) >0) {$btabpath = $opt_p;} else { $btabpath = $basedir; } # if no given for path to btab file, make default = base directory
if ((length($opt_t) >0) && (-e "$btabpath/$opt_t" == 1) && (-z "$btabpath/$opt_t" == 0)) {$btabfile = "$basedir/$opt_t";} else { &option_help; } # if no value for option t (name of btab file), quit with help menu
if ((length($opt_i) >0)  && (-e "$basedir/$opt_i" == 1) && (-z "$basedir/$opt_i" == 0)) {$infofile = "$basedir/$opt_i";} else { &option_help; } # if no value for option t (name of info file), quit with help menu
if ((length($opt_m) >0)  && (-e "$basedir/$opt_m" == 1) && (-z "$basedir/$opt_m" == 0)) {$hmmfile = "$basedir/$opt_m";} else { $hmmfile = undef;} # if no value for option m (name of hmm file), quit with help menu
if (length($opt_w) >0) {$window = $opt_w;} else { $window = 10000; } # if no window size given, make default "10000"
if (length($opt_s) >0) {$step = $opt_s;} else { $step = 5000; } # if no step size given, make default "5000"
if ($opt_E) {$evalue = $opt_E;} else { $evalue = 0.000001; } # if no E-value given, make default 0.000001
if ($opt_H) {$default_hpw = $opt_H;} else {$default_hpw = 4;} # if the # of hits per window is not defined, default = 4
if ($opt_a) {$defined_asmbl_id = $opt_a;} else { $defined_asmbl_id = undef; } # if assemble_id is not defined, then make $defined_asmbl_id undefined
if ((length($opt_r) >0) && (-e "$basedir/$opt_r" == 1) && (-z "$basedir/$opt_r" == 0)) {$tRNA_data = "$basedir/$opt_r";} else {$tRNA_data = undef;} # if file_name for tRNAscan-SE given, use it, if not, make $tRNA_data undefined
if ((length($opt_n) >0) && (-e "$basedir/$opt_n" == 1) && (-z "$basedir/$opt_n" == 0)) {$tmRNA_data = "$basedir/$opt_n";} else {$tmRNA_data = undef;} # if file_name for Aragorn data given, use it, if not, make $tmRNA_data undefined
if (length($opt_B) >0) {$asmbly_file_path = $opt_B;} else { $asmbly_file_path = $basedir; }
if ((length($opt_A) >0) && (-e "$asmbly_file_path/$opt_A" == 1) && (-z "$asmbly_file_path/$opt_A" == 0)) {$asmbly_file = "$asmbly_file_path/$opt_A"; $project_prefix = $opt_A; $project_prefix =~ s/\..*$//;} else {$asmbly_file = undef; $project_prefix = undef;} # if file_name for .1con given, use it, if not, make $file_name undefined (defined later)
if (($opt_F eq "M") || ($opt_F eq "m"))  {$search_method = "mummer";}
elsif (($opt_F eq "F") || ($opt_F eq "f"))  {$search_method = "fasta33";}
elsif (($opt_F eq "B") || ($opt_F eq "b"))  {$search_method = "blast";}
else { $search_method = "blast"; }
if ($opt_S) {$strict = 1;} else {$strict = 0;} # if strict mode (print only regions with core HMM hits or Mu-like and > 10 Kbp)
if ($opt_d) {$DEBUG = 1;} else {$DEBUG = 0;} # Debug mode.

############## Declare variables #################
my %fighash = (); # gets cleared with each assembly
my %hithash = (); # 1D, key = feat_name
my %rehash = (); # 2D, key1 = asmbl_id, key2 = end5
my %searchhash = (); # 2D, key1 = asmbl_id, key2 = end5
my %phage_hash = (); # gets cleared with each assembly
my %exclude_hash = (); # list of accessions to exclude
my %ok_comnames = ();
my %DB_info = (); # key is the BLAST DB tag, values are full name of phage and the taxonomy
my %exclude_hash = (); # list of phage accessions to NOT include
my %asmbl_idhash = ();
my %assembly_hash = (); # gets cleared with each assembly
my %tRNA_hash = ();
my %serine_HMM_hash = ('PF00239' => 1);
my %tyrosine_HMM_hash = ('PF00589' => 1,
                         'PF02899' => 1);
my %lists_hash = (); # stores the accessions for head morphogenesis Small and Large Terminases, Portal, Protease, Scaffolding, and Major Capsid proteins
my %HMMs_hash = (); # stores the curated HMM accessions
my $home = $ENV{'HOME'};
my $write_dir = "";
my @genomearray = ();
my $HMMversion = ""; # new variable to store the HMM version (currently 2 or 3)
my @hmm_data = ();
my $logfile = "$basedir/phage_phinder";
my $phome = "$home/phage_finder_v$version";
my $comfile = "$phome/phage_com_names_combo.txt";
my @filenames = ();
my $DBfile = "$phome/DB/phage_10_02_07_release.crib";
my $num_contigs = "";
my $phage_regions = "";
my $asmbly_status = 0; # boolean 0 = undefined, missing or inconsistant assembly data, 1 = ok assembly data
my $prefix = "";
my $hitsperwindow = $default_hpw; # initialize the hitsperwindow variable

sub find_regions {

  my ($asmbl_id,$hitsperwindow,$figref,$hitref,$phageref,$DEBUG) = @_;
  my $regions = 0;
  my $max_peak = 0;
  my $newcluster = 0;
  my $end5 = "";
  my $fig = "";
  my $hold_fig = "";

  foreach $fig (sort {$a <=> $b} keys %{$figref}) {
      if (($newcluster == 1) && ($figref->{$fig}->{'counts'} < ($hitsperwindow-1)))  { # if really the end of the region and not a window with 0 hits in the middle 
         # note this > 35000 condition is a temporary fix.  Really want both peaks processed and then a decision made at the end to determine if 1 or 2 phages (piggy back or same phage)
	 $regions++ if ((!defined($phageref->{$regions}->{'seed_end5'})) || (($end5 - $phageref->{$regions}->{'seed_end5'}) > 35000)); # we don't have 2 peaks within 35000 bp
         $phageref->{$regions}->{'seed_end5'} = $end5;
         $figref->{$hold_fig}->{'peak_number'} = $regions;
         $figref->{$hold_fig}->{'peak_value'} = $max_peak;
         $max_peak = 0; # reset max peak indicator to zero
         $newcluster = 0;
         print "region = $regions, end5 = $end5, seed_end5 = $phageref->{$regions}->{'seed_end5'}\n" if ($DEBUG);
      }
      else {
	  $newcluster = 0;
      }
      if (($figref->{$fig}->{'counts'} >= $hitsperwindow) && ($figref->{$fig}->{'counts'} > $max_peak))  {  #### Default = 4, ver 1.9beta changed to user defined ?
	  $max_peak = $figref->{$fig}->{'counts'};
	  $end5 = $hitref->{${$figref->{$fig}->{'featnames'}}[0]}->{'end5'};
          print "HITSperWINDOW = $hitsperwindow\n" if ($DEBUG);
	  print "FOUND asmbl_id ($asmbl_id): $regions - max_peak = $max_peak, and end5 = $end5\n" if ($DEBUG);
      }
      elsif (($max_peak > 0) && ($figref->{$fig}->{'counts'} == 0)) {# && ($figref->{($fig+$step)}->{'counts'} == 0))  {
	  $newcluster = 1;
	  $hold_fig = $fig;
      }
  }
  if (($max_peak > 0) && ($newcluster == 0))  {  # case where we run out of sequence and never got 0 phage hits (found bug with scaffold searches)
    $regions++;
    $phageref->{$regions}->{'seed_end5'} = $end5;
    print "FIX:  $end5\t$phageref->{$regions}->{'seed_end5'}\t$regions\n" if ($DEBUG);
  }
  delete $phageref->{"0"}; # delete the zero key
  return ($regions);
}

sub find_5prime_end {

    my ($n,$asmbl_id,$aref,$phageref,$hitref,$tRNAref,$reref,$okref,$figref,$DEBUG) = @_;
    my $i = "";
    my $lastkey = $aref->[$#{$aref}]+1; # start $lastkey as the 3' most coordinate
    my $hold = "";
    my $j = "";
    my $find = "";
    my $rRNAdist = 2000;
    my $feat_name = "";
    
    print "find_5prime_end: #phage_regions = $phage_regions\n" if ($DEBUG);
    for ($i = $#{$aref}; $i >= 0; $i--)  {
        print "|====>i = [$i]<====|\n" if ($DEBUG);
	if (($aref->[$i] <= $phageref->{$n}->{'seed_end5'}) && ($n >= 1) && ($aref->[$i] < $lastkey))  {
            print "------->n = $n,i = $i, aref = $aref->[$i]<-------\n" if ($DEBUG);
	    my $lm = nlowmult($window, $aref->[$i])-$window;
	    print "find_5prime_end: n = $n\tworking on $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}, clean_name = <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}>, lowmult = \[$lm\]...\n" if ($DEBUG);
	    print "COUNTS: $figref->{nlowmult($window, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]\n" if ($DEBUG);
	    print "nextCOUNTS: $figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'}\tnextNAME: $figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]\n" if ($DEBUG);
            if (exists $hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hmm'}) {
		print "find_5prime_end: HMM5(1a): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
                $phageref->{$n}->{'last5ORF'} = $lastkey;
	    }
	    # if this ORF has a hit in the Phage db, then include
	    elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hit'} == 1) {  
		print "find_5prime_end: HIT5(1a): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
                $phageref->{$n}->{'last5ORF'} = $lastkey;
	    }
	    # if tRNA or tmRNA is present and the distance between last orf and tRNA is < 2 Kbp (to avoid getting those between rRNAs), include $phageref->{$key}->{'tRNA'}{$t}
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i]}) && (($lastkey - $aref->[$i]) < $rRNAdist))  { 
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i]} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
		print "find_5prime_end: tRNA5: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-1]}) && (($lastkey - $aref->[$i-1]) < $rRNAdist) && ($i-1 >= 0))  { 
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i-1]} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
		$lastkey = $aref->[$i-1];
		print "find_5prime_end: tRNA5bypass(1a): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-1a): $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-2]}) && (($lastkey - $aref->[$i-2]) < $rRNAdist) && ($i-2 >= 0) ) { 
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i-2]} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
		$lastkey = $aref->[$i-2];
		print "find_5prime_end: tRNA5bypass(1b): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-1b): $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-2b): $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'}\n" if ($DEBUG);
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i-3]}) && (($lastkey - $aref->[$i-3]) < $rRNAdist) && ($i-3 >= 0)) { 
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i-3]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i-3]} = $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'};
		$lastkey = $aref->[$i-3];
		print "find_5prime_end: tRNA5bypass(1c): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-1c): $reref->{$asmbl_id}{$aref->[$i-1]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-2c): $reref->{$asmbl_id}{$aref->[$i-2]}->{'featname'}\n" if ($DEBUG);
		print "find_5prime_end: tRNA5bypass(-3c): $reref->{$asmbl_id}{$aref->[$i-3]}->{'featname'}\n" if ($DEBUG);
	    }
	    # if ok common name and there are still hits in the current window, then include
	    elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
		  ($figref->{nlowmult($window, $aref->[$i])}->{'counts'} >= 1) && 
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
		{  
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last5ORF'} = $lastkey;
                    print ">>$figref->{nlowmult($window, $aref->[$i])}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]<<\n" if ($DEBUG);
		    print "find_5prime_end: OKcom_name5: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		}

	    # if ok common name and there are still hits in the next window, then include
            ### special note:  changed counts from >= 2 to >= 1 for the Campylobacter jejuni 84-25 genome since there were 2 peaks separated with a single hit in one window
	    elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
		  ($figref->{nlowmult($window, $aref->[$i])-$window}->{'counts'} >= 1) && 
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])-$window}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
		{  
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last5ORF'} = $lastkey;
                    print ">>$figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]<<\n" if ($DEBUG);
		    print "find_5prime_end: OKcom_name5c: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		}

	    # case where >= 3 hits in next window, but no clean name (bad annotation perhaps or something novel

	    elsif (($figref->{nlowmult($window, $aref->[$i])-$step}->{'counts'} >= 3) && 
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])-$step}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))
		{  
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last5ORF'} = $lastkey;
                    print ">>$figref->{nlowmult($window, $aref->[$i])}->{'counts'}\t$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]<<\n" if ($DEBUG);
		    print "find_5prime_end: OKcom_name5b: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		}

#        # if there are at least 2 hits in the current step and end5 of first featname having hit in current window is less than end 5 of current orf, include
	    elsif (($figref->{nlowmult($step, $aref->[$i])}->{'counts'} >= 2) && 
		   ($hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]}->{'end5'} < $aref->[$i]))  {   
		$hold = $hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]}->{'end5'}; # store coordinate of 5' most hit in current window
		$phageref->{$n}->{'memberhash'}{$hold}->{'featname'} = $reref->{$asmbl_id}{$hold}->{'featname'};
		$j = $i;
		until ($hold == $aref->[$j])  { # fill in the other featnames until hit reached          
		    $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
		    print "find_5prime_end: $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'}\n" if ($DEBUG);
		    if ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j-1] - $aref->[$j]) < $rRNAdist)) { 
			$phageref->{$n}->{'tRNA'}{$aref->[$j]} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
		    }
		    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j-1] - $aref->[$j]) >= $rRNAdist)) { # likely within a rRNA operon, remove
			$hold = $lastkey;
			until ($hold == $aref->[$j])  {
			    delete $phageref->{$n}->{'memberhash'}{$aref->[$j]};
                            delete $phageref->{$n}->{'tRNA'}{$aref->[$j]};
			    $j++;
			}
			    $feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
			if ($lastkey < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
                            if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 5' end
				$phageref->{$n}->{'5prime'} = $phageref->{$n}->{'last5ORF'};
			    }
                            else {
				$phageref->{$n}->{'5prime'} = $lastkey; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
			    }
			}
                        else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
			    $phageref->{$n}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
			}
			$phageref->{$n}->{'ORF5'} = $reref->{$asmbl_id}{$lastkey}->{'featname'}; # store feat_name of 5' most ORF
			print "find_5prime_end_1: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}\n" if ($DEBUG);
			$n = $n - 1; # decrement phage region key counter when finished with phage region
			if ($n == 0)  {return;}
                        next;
		    }
		    $j--;
		}
		$lastkey = $hold;
                $phageref->{$n}->{'last5ORF'} = $lastkey if (!exists $phageref->{$n}->{'tRNA'}{$lastkey});
	    }

	    # specifically check for the word "transposase" in the com_name - these are ok and are often not the first word
            elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'com_name'} =~ /transposase/) {
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
                $phageref->{$n}->{'last5ORF'} = $lastkey;
		print "find_5prime_end: OKtransposase5: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
	    }

	    else {
                ###### added 10/14/05 to include a tRNA that was target when there are no db hits leading to it up to 10 kb
                my $pos = $i;
                my $end5 = $aref->[$pos];
                ## added ($pos < 0) condition 05/14/07 because if phage region was at 5prime end of a contig, this would loop indefinately
                until ((($lastkey-$end5) > 10000) || ($pos < 0)){ # only check within 10 kb of lastkey
		    $pos--;
                    $end5 = $aref->[$pos];
                    if (exists $tRNAref->{$asmbl_id}{$end5}) {
			$phageref->{$n}->{'memberhash'}{$end5}->{'featname'} = $reref->{$asmbl_id}{$end5}->{'featname'};
			$phageref->{$n}->{'tRNA'}{$end5} = $reref->{$asmbl_id}{$end5}->{'featname'};
			$lastkey = $end5;
			print "find_5prime_end: tRNAextend: $reref->{$asmbl_id}{$end5}->{'featname'}\n" if ($DEBUG);
                    }
                }		    
		$feat_name = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
		if ($lastkey < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
		    $phageref->{$n}->{'5prime'} = $phageref->{$n}->{'last5ORF'}; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
                    print "find_5prime_end: -----> feat_name = $feat_name, last5ORF = $phageref->{$n}->{'last5ORF'}, 5prime end of region = $phageref->{$n}->{'5prime'}\n" if ($DEBUG);
		}
		else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
		    $phageref->{$n}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
		    print "find_5prime_end: <----- feat_name = $feat_name, last5ORF = $phageref->{$n}->{'last5ORF'}, 5prime end of region = $phageref->{$n}->{'5prime'}\n" if ($DEBUG);
	        }
		$phageref->{$n}->{'ORF5'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
		print "find_5prime_end_2: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}\n" if ($DEBUG);
		$n = $n - 1; # decrement phage region key counter when finished with phage region
		if ($n == 0)  {return;}
	    }
	    
	}
    }
    if ($n == 1)  {  # weird case when phage is at the end of a molecule (when plasmid in genome project is actually the replicative form of an induced prophage
	$feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	if ($phageref->{$n}->{'last5ORF'} < $hitref->{$feat_name}->{'end3'})  {  # if 5 prime ORF ----->
	    $phageref->{$n}->{'5prime'} = $phageref->{$n}->{'last5ORF'}; # if find "housekeeping gene" quit the search and record the 5' end as the last good orf
	}
	else {   # <---- make 3' end of 5' most ORF the 5' boundary of phage region
	    $phageref->{$n}->{'5prime'} = $hitref->{$feat_name}->{'end3'};
	}
	$phageref->{$n}->{'ORF5'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last5ORF'}}->{'featname'};
	print "find_5prime_end_3: ==> Region $n 5prime end = $phageref->{$n}->{'5prime'}\tORF = $phageref->{$n}->{'ORF5'}\n" if ($DEBUG);
	$n = $n - 1; # decrement phage region key counter when finished with phage region
    }
}


sub find_3prime_end {

    my ($max,$asmbl_id,$aref,$phageref,$hitref,$tRNAref,$reref,$okref,$figref,$DEBUG) = @_;
    my $i = "";
    my $j = "";
    my $n = 1;
    my $hold = "";
    my $lastkey = "";
    my $find = "";
    my $rRNAdist = 2000;
    my $feat_name = "";
    
    for $i ( 0 .. $#{$aref} )  {
	if (($aref->[$i] >= $phageref->{$n}->{'seed_end5'}) && ($n <= $max)  && ($aref->[$i] > $lastkey))  {
	    print "find_3prime_end: n = $n\n" if ($DEBUG);  
            print "CURRENT FEAT_NAME:  $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}, clean_name = <$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}>\n" if ($DEBUG);
	    print "COUNTS: $figref->{nlowmult($window, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[0]..$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])+$step}->{'featnames'} }]\n" if ($DEBUG);
	    print "currentSTEPcounts:  $figref->{nlowmult($step, $aref->[$i])}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]\n" if ($DEBUG);
	    print "nextSTEPcounts:  $figref->{nlowmult($step, $aref->[$i]+$step)}->{'counts'}\tNAME: $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'}[0]..$figref->{nlowmult($step, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i]+$step)}->{'featnames'} }]\n" if ($DEBUG);
            print "nextWINDOWcounts: $figref->{nlowmult($window, $aref->[$i])+$window}->{'counts'}\tNAME: $figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'}[0]..$figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i]+$window)}->{'featnames'} }]\n" if ($DEBUG);
	    if (exists $hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hmm'}) {
		print "find_3prime_end:HMM3: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
                $phageref->{$n}->{'last3ORF'} = $lastkey;
	    }  
	    # if this ORF has a hit in the Phage db, then include
	    elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'hit'} == 1) {
		print "find_3prime_end: HIT3(1a): $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
                $phageref->{$n}->{'last3ORF'} = $lastkey;
	    }
	    # if tRNA or tmRNA is present, include
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i]}) && (($aref->[$i] - $lastkey) < $rRNAdist)) { 
		print "find_3prime_end: tRNA3: $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i]} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+1]}) && (($aref->[$i+1] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
		print "find_3prime_end: tRNA3bypass(1a): $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+1a): $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i+1]} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
		$lastkey = $aref->[$i+1];
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+2]}) && (($aref->[$i+2] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
		print "find_3prime_end: tRNA3bypass(1b): $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+1b): $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+2b): $aref->[$i+2]\t$reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i+2]} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
		$lastkey = $aref->[$i+2];
	    }
	    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$i+3]}) && (($aref->[$i+3] - $lastkey) < $rRNAdist) && ($i+1 <= $#{$aref})) {
		print "find_3prime_end: tRNA3bypass(1c): $aref->[$i]\t$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+1c): $aref->[$i+1]\t$reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+2c): $aref->[$i+2]\t$reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'}\n" if ($DEBUG);
		print "find_3prime_end: tRNA3bypass(+3c): $aref->[$i+3]\t$reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+1]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+1]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+2]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+2]}->{'featname'};
		$phageref->{$n}->{'memberhash'}{$aref->[$i+3]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'};
		$phageref->{$n}->{'tRNA'}{$aref->[$i+3]} = $reref->{$asmbl_id}{$aref->[$i+3]}->{'featname'};
		$lastkey = $aref->[$i+3];
	    }
            # if ok common name and there are still hits in the current window, then include
	    elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
		  (($figref->{nlowmult($window, $aref->[$i])}->{'counts'}) >= 1) &&
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
		{ 
		    print "find_3prime_end: OKcom_name3: ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}) $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG); 
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last3ORF'} = $lastkey;
		}
            
            # if ok common name and there are still hits in the next window, then include
	    elsif (($okref->{$hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'clean_name'}} == 1) &&
		  (($figref->{nlowmult($window, $aref->[$i])+$window}->{'counts'}) >= 1) &&
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])+$window}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
		{ 
		    print "find_3prime_end: OKcom_name3b: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG); 
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last3ORF'} = $lastkey;
		}

	    # if not good com_name and there are still >=2 (changed from 3) hits in the next step, then include
	    elsif ((($figref->{nlowmult($window, $aref->[$i])+$step}->{'counts'}) >= 2) &&
		  ($hitref->{$figref->{nlowmult($window, $aref->[$i])+$step}->{'featnames'}[$#{ $figref->{nlowmult($window, $aref->[$i])+$step}->{'featnames'} }]}->{'end5'} > $aref->[$i]))
		{ 
		    print "find_3prime_end: OKcom_name3b: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG); 
		    $phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		    $lastkey = $aref->[$i];
                    $phageref->{$n}->{'last3ORF'} = $lastkey;
		}
	    
	    # if there is a hit in the current step and end3 of the last featname having hit in current window is less than end 5 of current orf, include
	    elsif ((($figref->{nlowmult($step, $aref->[$i])}->{'counts'}) >= 2) && 
		   ($hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]}->{'end5'} > $aref->[$i])) {  
		$hold = $hitref->{$figref->{nlowmult($step, $aref->[$i])}->{'featnames'}[$#{ $figref->{nlowmult($step, $aref->[$i])}->{'featnames'} }]}->{'end5'};
                # if there are more than 3 ORFs until the next hit, then bail
		$phageref->{$n}->{'memberhash'}{$hold}->{'featname'} = $reref->{$asmbl_id}{$hold}->{'featname'};
		$j = $i;
		until ($hold == $aref->[$j])  { 
                    print "find_3prime_end:  fill to hit in next step [ $reref->{$asmbl_id}{$aref->[$j]}->{'featname'} ]\n" if ($DEBUG);          
		    $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
		    if ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j] - $aref->[$j-1]) < $rRNAdist)) { 
			$phageref->{$n}->{'tRNA'}{$aref->[$j]} = $reref->{$asmbl_id}{$aref->[$j]}->{'featname'};
			print "find_3prime_end: STEP: $phageref->{$n}->{'memberhash'}{$aref->[$j]}->{'featname'}\n" if ($DEBUG);
		    }
		    elsif ((exists $tRNAref->{$asmbl_id}{$aref->[$j]}) && (($aref->[$j] - $aref->[$j-1]) >= $rRNAdist)) { # likely within a rRNA operon, remove
			$hold = $lastkey;
			until ($hold == $aref->[$j])  {
			    delete $phageref->{$n}->{'memberhash'}{$aref->[$j]};
                            delete $phageref->{$n}->{'tRNA'}{$aref->[$j]};
			    $j--;
			}
			$feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
			if ($lastkey > $hitref->{$feat_name}->{'end3'})  {  # if 3 prime ORF <-----
                            if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 3' end
                                $phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'};
			    }
			    else {
				$phageref->{$n}->{'3prime'} = $lastkey; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
			    }
			}
			else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
			    $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
			}
			$phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$lastkey}->{'featname'};
			print "find_3prime_end: STEP ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
			$n = $n + 1; # increment phage region key counter when finished with phage region
			if ($n == $max+1) {return;}
                        next;
		    }
		    $j++;
		}
		$lastkey = $hold;
		$phageref->{$n}->{'last3ORF'} = $lastkey if (!exists $phageref->{$n}->{'tRNA'}{$lastkey});
	    }	    
	    # specifically check for the word "transposase" in the com_name - these are ok and are often not the first word
            elsif ($hitref->{$reref->{$asmbl_id}{$aref->[$i]}->{'featname'}}->{'com_name'} =~ /transposase/) {
		print "find_3prime_end: OKtransposase3: $reref->{$asmbl_id}{$aref->[$i]}->{'featname'}\n" if ($DEBUG);
		$phageref->{$n}->{'memberhash'}{$aref->[$i]}->{'featname'} = $reref->{$asmbl_id}{$aref->[$i]}->{'featname'};
		$lastkey = $aref->[$i];
	    }
	    else {
########### added 05/14/07 to include a tRNA that was target when there are no db hits leading to it up to 10 kb
                my $pos = $i;
                my $end5 = $aref->[$pos];
                ## added ($pos > $#{$aref}) condition 05/14/07 because if phage region was at 3prime end of a contig, this would loop indefinately
                until ((($end5-$lastkey) > 10000) || ($pos > $#{$aref})){ # only check within 10 kb of lastkey
		    $pos++;
                    $end5 = $aref->[$pos];
                    if (exists $tRNAref->{$asmbl_id}{$end5}) {
			$phageref->{$n}->{'memberhash'}{$end5}->{'featname'} = $reref->{$asmbl_id}{$end5}->{'featname'};
			$phageref->{$n}->{'tRNA'}{$end5} = $reref->{$asmbl_id}{$end5}->{'featname'};
			$lastkey = $end5;
			print "find_3prime_end: tRNAextend: $reref->{$asmbl_id}{$end5}->{'featname'}\n" if ($DEBUG);
                    }
                }		    
###########
		$feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
		if ($lastkey > $hitref->{$feat_name}->{'end3'})  {  # if 3 prime ORF <-----
		    if (exists $phageref->{$n}->{'tRNA'}{$lastkey}) { # case where tRNA is at 3' end
			$phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'};
		    }
		    else {
			$phageref->{$n}->{'3prime'} = $lastkey; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
		    }
		}
		else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
		    $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
		}
		$phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$phageref->{$n}->{'last3ORF'}}->{'featname'};
		print "find_3prime_end: ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
		$n = $n + 1; # increment phage region key counter when finished with phage region
		if ($n == $max+1) {return;}
	    }  
	}
    }
# weird case when phage is at the end of a molecule (when plasmid in genome project is actually the replicative form of an induced prophage
    if ($n <= $max)  {  
	$feat_name = $reref->{$asmbl_id}{$lastkey}->{'featname'};
	if ($phageref->{$n}->{'last3ORF'} == $phageref->{$n}->{'5prime'})  { 
	    $phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
	    $phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$phageref->{$n}->{'3prime'}}->{'featname'};
	}
	# if find "housekeeping gene" quit the search and record the 3' end as the last good orf
	else {
	    if ($phageref->{$n}->{'last3ORF'} > $hitref->{$feat_name}->{'end3'})  {  # if 3 prime ORF <-----
		$phageref->{$n}->{'3prime'} = $phageref->{$n}->{'last3ORF'}; # if find "housekeeping gene" quit the search and record the 3' end as the last good orf
	    }
	    else {   # ----> make 3' end of 5' most ORF the 3' boundary of phage region
		$phageref->{$n}->{'3prime'} = $hitref->{$feat_name}->{'end3'};
	    }
	    $phageref->{$n}->{'ORF3'} = $reref->{$asmbl_id}{$phageref->{$n}->{'3prime'}}->{'featname'};
	}
	print "find_3prime_end: ==> Region $n 3prime end = $phageref->{$n}->{'3prime'}\n" if ($DEBUG);
	$n = $n + 1; # increment phage region key counter when finished with phage region
    }
}

sub find_att_sites { # try to determine the attachment site for each phage region

    my ($prefix,$asmbl_id,$search_method,$ahref,$phageref,$hitref,$reref,$tRNAref,$aref,$figref,$window,$step,$hitsperwindow,$okref,$DEBUG) = @_;
    my @k = ();
    my @grase = ();
    my $key = ""; # stores keys of %phage_hash (phage number)
    my $yek = ""; # stores keys of %memberhash (feat_names)
    my $i = ""; # reused loop var
    my $phage_5prime = "" ; # predicted beginning of phage region
    my $phage_3prime = ""; # predicted end of phage region
    my $halfsize = ""; # middle of phage region coordinate
    my $Bint = ""; # leftmost end of left most integrase (either 5' or 3' of orf, but 5' relative to phage region)
    my $Eint = ""; # rightmost end of right most integrase (either 5' or 3' of orf, but 3' relative to phage region)
    my $end5 = ""; # 5' end of integrase
    my $end3 = ""; # 3' end of integrase
    my $head = ""; # 400 bp region that is outward-facing from Bint or Eint
    my $tail = ""; # region that is the remainder of the phage region plus 5000 bp for slop
    my $size = ""; # size from Bint or Eint until the end of the predicted phage region 
    my $size_head = ""; # how large of head to pull
    my $size_tail = "";  # how large of tail to pull
    my $start_head = ""; # beginning of head region to pull
    my $start_tail = ""; # beginning of tail region to pull
    my $fivetend5 = "";
    my $fivetend3 = "";
    my $threetend5 = "";
    my $threetend3 = "";
    my $tRNA5 = "";
    my $tRNA3 = "";
    my $t = "";
    my $found_att = 0;
    my $int = 0;
    my $feat_name = "";
    my $hold_attB = ""; # placeholder for first attB found (beginning of region with tRNA att site
    my $hold_attE = ""; # placeholder for first attE found (end of region with tRNA att site
    my $start = ""; # contains coordinate in genomearray to start checking for target site
    my $finish = ""; # contains coordinate in genomearray to end the checking for target site
    my $lenB = ""; # stores the length of the att site for beginning int
    my $lenE = ""; # stores the length of the att site for ending int
    my $B_region_length = "";
    my $E_region_length = "";
    my $Bextension = 15000; # set size of extension (additional sequence past initial boundaries to search for TSD)
    my $extension = "";

###########################
    local *beginning_sub = sub {  # create a function local to sub find_att_sites (a nested subroutine) to minimize duplication of identical code and reusing vars
	($key) = @_;  # passed the $key variable because the local function was unable to see it otherwise
        $extension = $Bextension;
	$size = abs($Bint - $phage_3prime);
	$size_head = 400;
	$start_head = $Bint-$size_head;
        if ($start_head < 0) { # when working with contigs, we don't want to be running off the 5' end of the contig [ 05/15/07 ]
	    $start_head = 1;
            $size_head = $Bint;
	}
	$start_tail = $phage_3prime-(round($size*0.2));
	$size_tail = $extension+(round($size*0.2));
        if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # added 05/15/07 to make sure we don't run off the 3' end of a contig
	    $size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
            $extension = 0;
	}
	print "find_att_sites:beginning_sub: phage 5 prime end = $phage_5prime\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: phage 3 prime end = $phage_3prime\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size = $size\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: start_head = $start_head\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size_head = $size_head\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: start_tail = $start_tail\n" if ($DEBUG);
	print "find_att_sites:beginning_sub: size_tail = $size_tail\n" if ($DEBUG);
	($found_att, $hold_attB, $hold_attE) = &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"+","phage",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize);
	if ($found_att)  { # calculate the distance from outermost integrase coord to outermost att coord
	    $phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
	}
    };
############################
    local *ending_sub = sub {  # create a function local to sub find_att_sites (a nested subroutine) to minimize duplication of identical code and reusing vars
	($key) = @_; # passed the $key variable because the local function was unable to see it otherwise
	$extension = $Bextension;
	$size = abs($Eint - $phage_5prime);
	$size_head = 400;
        $start_head = $Eint;
        if (($start_head+$size_head) > $ahref->{$asmbl_id}->{'length'}) {  #05/15/07 to make sure we don't run off the 3' end of a contig/scaffold
	    $size_head = $ahref->{$asmbl_id}->{'length'}-$start_head;
        }
	$start_tail = $Eint-1-($size+$extension);
        if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
	    $start_tail = 1;
            $extension = 0;
	}
        $size_tail = $extension+(round($size*0.2));
	print "find_att_sites:ending_sub: phage 5 prime end = $phage_5prime\n" if ($DEBUG);
        print "find_att_sites:ending_sub: extension = $extension\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size = $size\n" if ($DEBUG);
	print "find_att_sites:ending_sub: start_head = $start_head\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size_head = $size_head\n" if ($DEBUG);
	print "find_att_sites:ending_sub: start_tail = $start_tail\n" if ($DEBUG);
	print "find_att_sites:ending_sub: size_tail = $size_tail\n" if ($DEBUG);
	($found_att, $hold_attB, $hold_attE) = &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"-","phage",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize);
	if ($found_att)  { # calculate the distance from outermost integrase coord to outermost att coord
	    $phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
	}
    };
############################
    local *tRNA_sub = sub {  # create a function local to sub find_att_sites to find tRNA TSDs
	($key,$Bint,$Eint) = @_; # passed the $key variable because the local function was unable to see it otherwise
        $extension = $Bextension;
	$fivetend5 = $tRNA5;
	$fivetend3 = $tRNAref->{$asmbl_id}{$tRNA5}->{'end3'};
	$threetend5 = $tRNA3;
	$threetend3 = $tRNAref->{$asmbl_id}{$tRNA3}->{'end3'};
	print "find_att_sites:tRNA_sub: fivetend5: $fivetend5\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: fivetend3: $fivetend3\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: threetend5: $threetend5\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: threetend3: $threetend3\n" if ($DEBUG);
	print "find_att_sites:tRNA_sub: halfsize: $halfsize\n" if ($DEBUG);
	# only viably possiblities are -----> at beginning or <----- at end of region
	if (defined $fivetend5)  { # if tRNA is on the 5' side of the region
	    if ($fivetend5 < $fivetend3) { # direction ---> and on 5' end of phage region
		$start_head = $fivetend5-1;
                if ($start_head < 0) { # 05/15/07
		    $start_head = 1;
		}
		$start_tail = $fivetend3 + 1;
		$size_head = ($fivetend3 - $fivetend5) + 1; # pull complete tRNA and search this first
		$size_tail = ($phage_3prime - $start_tail) + 1 + $extension;
                if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # 05/15/07
		    $size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
		    $extension = 0;
		}
		print "find_att_sites:tRNA_sub: tRNA B----->\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		($found_att, $hold_attB, $hold_attE) = &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"+","phage\_$reref->{$asmbl_id}{$tRNA5}->{'featname'}",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize);
		if ($found_att == 1)  { # if match tRNA sequence, then check to see if we can expand on the match on 3' side
		    print "find_att_sites:tRNA_sub: second round expanded tRNA search...\n" if ($DEBUG);
		    $start_head = $fivetend5-1;
                    if ($start_head < 0) { # 05/15/07
			$start_head = 1;
		    }
		    $start_tail = $fivetend3 + 1 + 200;
		    $size_head = ($fivetend3 - $fivetend5) + 1 + 200; # pull complete tRNA + 200 extra nucleotides on 3' end
		    $size_tail = ($phage_3prime - $start_tail) + 1 + $extension + 200;
                    if (($start_tail+$size_tail) > $ahref->{$asmbl_id}->{'length'}) { # 05/15/07
			$size_tail = $ahref->{$asmbl_id}->{'length'}-$start_tail;
			$extension = 0;
		    }
		    print "find_att_sites:tRNA_sub: tRNA B----->\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		    &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"+","phage\_$reref->{$asmbl_id}{$tRNA5}->{'featname'}\_r2",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,$hold_attB,$hold_attE,$figref,$window,$step,$hitsperwindow,$halfsize);
		    $phageref->{$key}->{'target'} = $reref->{$asmbl_id}{$tRNA5}->{'featname'};
		}
		else { print "find_att_sites:tRNA_sub: NO att site found for $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);}
		}
	    else { print "find_att_sites:tRNA_sub: $reref->{$asmbl_id}{$tRNA5}->{'featname'} facing wrong direction <---- at beginning of phage region $key\n" if ($DEBUG);}
	    }
	if ((defined $threetend5) && ($found_att == 0)) { # if tRNA is on the 3' side of the region
	    if (($threetend5 > $threetend3) && ($found_att == 0) ) { # direction <---- and not 1 tRNA and 3' end of phage region and no att site found previously
		$start_head = $threetend3; # pull complete tRNA and search this first
		$size_head = ($threetend5 - $threetend3) + 1;
		$start_tail = $phage_5prime - $extension;
                if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
		    $start_tail = 1;
		    $extension = 0;
		}
		$size_tail = (($start_head-1) - $start_tail) + 1;
		print "find_att_sites:tRNA_sub: tRNA <-----E\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		($found_att, $hold_attB, $hold_attE) = &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"-","phage\_$reref->{$asmbl_id}{$tRNA3}->{'featname'}",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,"0","0",$figref,$window,$step,$hitsperwindow,$halfsize);
		if ($found_att == 1)  { # if match tRNA sequence, then check to see if we can expand on the match on 3' side
		    print "find_att_sites:tRNA_sub: second round expanded tRNA search...\n" if ($DEBUG);
		    $start_head = $threetend3-201; # pull 200 extra on 3' side to search
		    $size_head = ($threetend5 - $threetend3) + 1 + 200;
		    $start_tail = $phage_5prime - $extension;
		    if ($start_tail < 0) { # added 05/15/07 to make sure we don't run off the 5' end of a contig/scaffold
			$start_tail = 1;
			$extension = 0;
		    }
		    $size_tail = (($start_head-1) - $start_tail) + 1;
		    print "find_att_sites:tRNA_sub: tRNA <-----E\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_head = $start_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_head = $size_head\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: start_tail = $start_tail\n" if ($DEBUG);
		    print "find_att_sites:tRNA_sub: size_tail = $size_tail\n" if ($DEBUG);
		    &get_TSD($prefix,$asmbl_id,$hitref,$phageref,$reref,$aref,$DEBUG,"-","phage\_$reref->{$asmbl_id}{$tRNA3}->{'featname'}\_r2",$search_method,$key,$ahref,$start_head,$start_tail,$size_head,$size_tail,$hold_attB,$hold_attE,$figref,$window,$step,$hitsperwindow,$halfsize);
		    $phageref->{$key}->{'target'} = $reref->{$asmbl_id}{$tRNA3}->{'featname'};
		}
		else { print "find_att_sites:tRNA_sub: NO att site found for $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);}
		}
	    else { print "find_att_sites:tRNA_sub: $reref->{$asmbl_id}{$tRNA3}->{'featname'} facing wrong direction ----> at end of phage region $key\n" if ($DEBUG);}
	}
        if ($found_att)  {
	    if (($Bint > "") && ($Eint > "")) {
		if (($Bint - $phageref->{$key}->{'left'}) == ($phageref->{$key}->{'right'} - $Eint)) { # if for some reason there are 2 integrases and they have the same distance to the att sites
		    if (defined $fivetend5)  {
			$phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
		    }
		    else {
			$phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
		    }
		}
		elsif (($Bint - $phageref->{$key}->{'left'}) < ($phageref->{$key}->{'right'} - $Eint)) {
		    $phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
		}
		else  {
		    $phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
		}
	    }
	    elsif (($Bint > "") && ($Eint == "")) {
		$phageref->{$key}->{'att_distance'} = $Bint - $phageref->{$key}->{'left'};
	    }
	    elsif (($Bint == "") && ($Eint > "")) {
		$phageref->{$key}->{'att_distance'} = $phageref->{$key}->{'right'} - $Eint;
	    }
	}
    };
############################
    foreach $key (sort {$a <=> $b} keys %{$phageref}) {
	print "find_att_sites:MAIN LOOP: REGION = $key\n" if ($DEBUG);
	$phage_5prime = $phageref->{$key}->{'5prime'};
	$phage_3prime = $phageref->{$key}->{'3prime'};
        print "find_att_sites:  The 5prime end before att is $phage_5prime, the 3prime end is $phage_3prime\n" if ($DEBUG);
	$size = $phage_3prime - $phage_5prime;
	print "find_att_sites:MAIN LOOP: SIZE = $size\n" if ($DEBUG);
	$halfsize = ($size*.5) + $phage_5prime;
	print "find_att_sites:MAIN LOOP: HALFSIZE = $halfsize\n" if ($DEBUG);
	$tRNA5 = undef;
	$tRNA3 = undef;
	if ($phageref->{$key}->{'class'} eq "Mu-like")  { next; } # skip if the region is a Mu-like phage (no integrase and no tRNA specificity)
	print "find_att_sites:PREINNER LOOP: END5 = $end5, END3 = $end3\n" if ($DEBUG);
	if ($phageref->{$key}->{'int_HMM'} > 0)  { # of outermost integrases
	    
########## check for valid integrases (first and last ###############
	    $int = 1;

# new logic to find 5' and 3' integrases
            push (@grase, shift @{$phageref->{$key}->{'integrases'}}); # gets 5' most integrase
            push (@grase, pop @{$phageref->{$key}->{'integrases'}}); # gets 3' most integrase if present
            print "5GRASE>>>>>> $grase[0]\t$grase[1]\n" if ($DEBUG);
            if ($grase[0] eq $grase[1]) { pop @grase; } # if the same integrase listed, remove last occurance
            for $i (@grase) {
		$end5 = $hitref->{$i}->{'end5'};
                $end3 = $hitref->{$i}->{'end3'};
		print "find_att_sites:INNER LOOP: END5 = $end5, END3 = $end3\n" if ($DEBUG);
		if (($end5 < $end3) && ($end5 < $halfsize) && ($Bint == "")) { # int is 5'--->3' and is in 5' half of phage region  
		    $Bint = $end5; # define the 5' most part of the int gene
									     }
		elsif (($end5 > $end3) && ($end5 < $halfsize) && ($Bint == "")) { # int is 3'<---5' and is in 5' half of phage region
		    $Bint = $end3; # define the 5' most part of the int gene 
										}
		elsif (($end5 > $end3) && ($end5 > $halfsize)) { 
		    $Eint = $end5;
		}
		elsif (($end5 < $end3) && ($end5 > $halfsize)) {
		    $Eint = $end3;
		}                
	    }	    
	} 
        else  {
	    print "find_att_sites: NO integrase HMMs found\n" if ($DEBUG);
	}
	if ((exists $phageref->{$key}->{'tRNA'}) && ($int == 1))  { # tRNA is present and an integrase in the region, find tRNA TSD
	    if (keys %{$phageref->{$key}->{'tRNA'}} == 1) { # if only 1 key (tRNA)
		@k = keys %{$phageref->{$key}->{'tRNA'}};
		if ($k[0] < $halfsize) {
		    $tRNA5 = $k[0];
		    print "===> $k[0], $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);
		}
		else {
		    $tRNA3 = $k[0];
		    print "===> $k[0], $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);
		}
	    }
	    elsif (keys %{$phageref->{$key}->{'tRNA'}} > 1) { # multiple tRNAs in region 
		foreach $t (sort {$a <=> $b} keys %{$phageref->{$key}->{'tRNA'}}) { # get most 5' tRNA
		    print "t: $t\ttRNA: $phageref->{$key}->{'tRNA'}{$t}\n" if ($DEBUG);
		    last if (defined($tRNA5) && (($t-$tRNA5) > 1000) && ($tRNA5 < $tRNAref->{$asmbl_id}{$tRNA5}->{'end3'})); # stop looking if we have a 5' tRNA and the distance between 2 tRNAs is > 1000 bp and tRNA must be in proper orientation
		    $tRNA5 = $t if ($t < $halfsize);
		    print "tRNA5===> $tRNA5, $reref->{$asmbl_id}{$tRNA5}->{'featname'}\n" if ($DEBUG);
		}
		foreach $t (sort {$b <=> $a} keys %{$phageref->{$key}->{'tRNA'}}) { # get most 3' tRNA
		    print "t: $t\ttRNA: $phageref->{$key}->{'tRNA'}{$t}\n" if ($DEBUG);
		    last if (defined($tRNA3) && (($tRNA3-$t) > 1000) && ($tRNA3 > $tRNAref->{$asmbl_id}{$tRNA3}->{'end3'}));
		    $tRNA3 = $t if ($t > $halfsize);
		    print "tRNA3===> $tRNA3, $reref->{$asmbl_id}{$tRNA3}->{'featname'}\n" if ($DEBUG);
		}
	    }
	    &tRNA_sub($key,$Bint,$Eint);
	    print "find_att_sites: TRNA found_att -----> tRNA5: $tRNA5 ($reref->{$asmbl_id}{$tRNA5}->{'featname'})\ttRNA3: $tRNA3 ($reref->{$asmbl_id}{$tRNA3}->{'featname'})\n"if ($DEBUG);
            print "find_att_sites: Bint = $Bint, Eint = $Eint\n" if ($DEBUG);
            print "find_att_sites: FOUND_att = $found_att, int = $int\n" if ($DEBUG);
	}
	if (($found_att == 0) && ($int == 1))  { # if no tRNA att found and there are integrases in region, proceed
	    
########## if no tRNA found, check for TSDs ###############
	    
	    if (($Bint > "") && ($Eint > "")) { # If there are more than one putative integrases present, check both and pick best att (1) target, 2)longest length)
		print "find_att_sites: Phage region $key has at least 2 integrases, at $Bint and $Eint\n" if ($DEBUG);
		print "find_att_sites: Checking Beginning integrase (1 of 2) for phage region $key begins at Beginning $Bint A\n" if ($DEBUG);
		&beginning_sub($key);
		$lenB = length($hitref->{"attL_$asmbl_id\_$key"}->{'annotation'});
		if (!defined ($phageref->{$key}->{'target'})) {
		    if ($lenB > 0) { delete $hitref->{"attL_$asmbl_id\_$key"}; } # removing old data so new data can be stored
		    print "find_att_sites: Checking Ending integrase (2 of 2) for phage region $key begins at End $Eint B\n" if ($DEBUG);
		    delete $hitref->{"attL_$asmbl_id\_$key"}; # delete data from beginning_sub
		    $phageref->{$key}->{'target'} = "";
		    print "hold_attB = $hold_attB\thold_attE = $hold_attE\n" if ($DEBUG);
		    delete $phageref->{$key}->{'memberhash'}{$hold_attB};  # remove old attB information
		    delete $phageref->{$key}->{'memberhash'}{$hold_attE};  # remove old attE information
                    $B_region_length = abs($hold_attB-$hold_attE);
		    &ending_sub($key);
                    $lenE = length($hitref->{"attL_$asmbl_id\_$key"}->{'annotation'});
                    $E_region_length = abs($hold_attB-$hold_attE);
                    print "Blength = $B_region_length\tElength = $E_region_length\n" if ($DEBUG);
		    if (($B_region_length > $E_region_length) || ($lenB > $lenE))  { # if length of phage region with B int is > or if the length of the att site is >
			print "find_att_sites: Going with beginning integrase (lenB = $lenB, lenE = $lenE) when 2 present\n" if ($DEBUG);
			delete $hitref->{"attL_$asmbl_id\_$key"}; # delete data from ending_sub
			$phageref->{$key}->{'target'} = "";
                        print "hold_attB = $hold_attB\thold_attE = $hold_attE\n" if ($DEBUG);
                        delete $phageref->{$key}->{'memberhash'}{$hold_attB};  # remove old attB information
			delete $phageref->{$key}->{'memberhash'}{$hold_attE};  # remove old attE information
			&beginning_sub($key);
		    }
		}
	    }
	    elsif (($Bint > "") && ($Eint == "")) { # If only one at the beginning, use it
		print "find_att_sites: The integrase to study for phage region $key begins at Beginning $Bint C\n" if ($DEBUG);
		&beginning_sub($key);
	    }
	    elsif (($Bint == "") && ($Eint > "")) { # If only one at the end, use it
		print "find_att_sites: The integrase to study for phage region $key begins at End $Eint D $Bint\n" if ($DEBUG);
		&ending_sub($key);
	    }
	    else { print "find_att_sites: No integrases for phage region $key\n" if (($DEBUG) && ($int == 0)); }
	}

	if ($found_att)  { # if an att site found and the new boundaries are beyond the original phage boundaries, then add ORFs
            print "FOUND_att!\n" if ($DEBUG);
	    if (($phageref->{$key}->{'left'} < $phageref->{$key}->{'5prime'}) || ($phageref->{$key}->{'right'} > $phageref->{$key}->{'3prime'}))  {
		print "EXTENDING REGION!\n" if ($DEBUG);
		my $rv = &add_ORFs_to_region($key,$asmbl_id,$phageref,$hitref,$reref,$aref,$okref,$DEBUG);
                if ($rv)  {
		    print "find_att_sites: House-keeping genes found, undoing att info...\n" if ($DEBUG);
		    $phageref->{$key}->{'5prime_att'} = "";
                    $phageref->{$key}->{'3prime_att'} = "";
                    $phageref->{$key}->{'left'} = "";
		    $phageref->{$key}->{'right'} = "";
                    $phageref->{$key}->{'att_distance'} = "";
                    $phageref->{$key}->{'target'} = "";
                    print "EXTENDING_ORFS: Direction = <$phageref->{$key}->{'direction'}>\n" if ($DEBUG);
		    $phageref->{$key}->{'direction'} = "+"; #reset to default (+)
		    if ($phageref->{$key}->{'5prime'} > $phageref->{$key}->{'3prime'}) {
			print "changing directions\n" if ($DEBUG);
			my $hold_5prime = $phageref->{$key}->{'5prime'};
			my $hold_3prime = $phageref->{$key}->{'3prime'};
			$phageref->{$key}->{'5prime'} = $hold_3prime;
			$phageref->{$key}->{'3prime'} = $hold_5prime;
		    }
                    $found_att = 0;
                    next;
		}
		
	    }
	    if ($phageref->{$key}->{'att_distance'} < 0)  { # if att_distance is negative, make zero (case where real int is not translated and found one outside of the region
                print "find_att_sites:  att_distance is $phageref->{$key}->{'att_distance'}!, so setting att_distance to zero\n" if ($DEBUG);
		$phageref->{$key}->{'att_distance'} = 0;
	    }
	}
	$found_att = 0; # reset so not found
	$int = 0;
	$Bint = "";
	$Eint = "";
	$end5 = "";
	$end3 = "";
        @grase = ();
    }
}

sub print_regions {

my ($write_dir,$asmbly_status,$asmbl_id,$hitsperwindow,$window,$step,$evalue,$max,$phageref,$aref,$ahref,$hitref,$reref,$dbref,$DEBUG) = @_; 
my $large = 0;
my $totalphagebp = 0;
my $percent = 0;
my $pico = 0;  # size < 2
my $nano = 0;  # size >= 2, < 5
my $micro = 0; # size >= 5, < 10
my $small = 0;
my $medium = 0;
my $size = 0;
my $size_att = 0;
my $sizer = 0;
my $totalpicobp = 0;
my $totalnanobp = 0;
my $totalmicrobp = 0;
my $totalsmallbp = 0;
my $totalmediumbp = 0;
my $percentsmall = 0;
my $percentmedium = 0;
my $hmm_hit_names = "";
my $end5 = "";
my $end3 = "";
my $featname = "";
my $printfeatname = "";
my $annotation = "";
my $header = "";
my $con_header = "";
my $region_header = "";
my $region2_header = "";
my $region_tag = "";
my $seq_tag = "";
my $key = ""; 
my $yek = ""; 
my $HMM_name = "";
my $small_open = 0;
my $med_open = 0;
my $large_open = 0;
my $phage_end5 = "";
my $phage_end3 = "";
my $phage_direction = "";
my $label = "";
my $bugtag = $ahref->{$asmbl_id}->{'genus'} . "_" . $ahref->{$asmbl_id}->{'species'} . "_" . $ahref->{$asmbl_id}->{'strain'} . "_" . $asmbl_id;
my $prefix = "$write_dir";
my $seqfileprefix = "$prefix/$bugtag";
my $pepfileprefix = "$prefix/$bugtag";
my $confileprefix = "$prefix/$bugtag";
my $print_files = undef;

local *print_sub = sub { # universal print subroutine within print_regions
  my ($key, $FH, $label, $count, $type, $phage_end5, $phage_end3, $phage_direction) = @_;
  my $hold_locus = "";
  my $print5 = 0;
  my $seq = "";
  my $pep = "";
  my $seq_header = "";
  my $i = "";
  my $B_pos = "";
  my $E_pos = "";
  my %Best_hits = (); # reset with each call to print_sub
  my $k = "";
  my $hit = "";
  my $best = "";
  my $shortag = "";
  my $totalhits = 0;
  $phageref->{$key}->{'ORFcnt'} = 0; # set the ORF count to zero

  if ($count > 1)  {
          print $FH "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
          select((select($FH), $= = 67, $- = 0)[0]);
        }
        else  {
          select((select($FH), $= = 67)[0]);
        }
  print "DIRECTION: <$phage_direction>\n" if ($DEBUG);
   foreach $k (sort {$a <=> $b} keys %{$phageref->{$key}->{'memberhash'}}) { # pre loop to get the Best Hit
       $featname = $phageref->{$key}->{'memberhash'}{$k}->{'featname'};
       print "BEST_HITS: pre_featname = $featname\n" if ($DEBUG);
       if ((exists $hitref->{$featname}->{'annotation'}) && ($featname !~ /^att/) && ($featname !~ /^tm?RNA/)){
	   if (($hitref->{$featname}->{'end5'} >= $phage_end5) && ($hitref->{$featname}->{'end5'} <= $phage_end3))  {
	       $Best_hits{$hitref->{$featname}->{'phage'}}->{'counts'}++;
               $totalhits++;
	       print "BEST_HITS:    post_featname: $featname\t$hitref->{$featname}->{'annotation'}\t$hitref->{$featname}->{'phage'}\n" if ($DEBUG);
	   }
       }
   }
   foreach $hit (sort {$Best_hits{$b}->{'counts'} <=> $Best_hits{$a}->{'counts'}} keys %Best_hits)  {
       $best = $hit;
       #print "$best ($Best_hits{$best}->{'counts'})\n";
       $shortag = $best;
       $shortag =~ s/^PHAGE_//;
       last;
   }
  if (length($phageref->{$key}->{'target'}) > 0) {
    print "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $phageref->{$key}->{'target'}. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
    $region_header = "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%, TARGET: $phageref->{$key}->{'target'}.";
##
    print TAB "$asmbl_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\t$size_att\t$label\t$type\t$hitref->{\"attR_$asmbl_id\_$key\"}->{'annotation'}\t$hitref->{\"attL_$asmbl_id\_$key\"}->{'annotation'}\t$phageref->{$key}->{'target'}\t$phageref->{$key}->{'gc'}\t$shortag\t";
##
  
  }
  elsif ($size_att > 1)  {
    print "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
    $region_header = "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} or ($phageref->{$key}->{'5prime_att'} to $phageref->{$key}->{'3prime_att'}) and is $size ($size_att) bp in size, gc%: $phageref->{$key}->{'gc'}%.";

##
    print TAB "$asmbl_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime_att'}\t$phageref->{$key}->{'3prime_att'}\t$size_att\t$label\t$type\t$hitref->{\"attR_$asmbl_id\_$key\"}->{'annotation'}\t$hitref->{\"attL_$asmbl_id\_$key\"}->{'annotation'}\tN.D.\t$phageref->{$key}->{'gc'}\t$shortag\t";
##

  }
  else  {
    print "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%. Best hit: $best ($Best_hits{$best}->{'counts'}/$totalhits)\n";
    $region_header = "$label $type region $key is from $phageref->{$key}->{'5prime'} to $phageref->{$key}->{'3prime'} and is $size bp in size, gc%: $phageref->{$key}->{'gc'}%.";

##
    print TAB "$asmbl_id\t$asmbl_idhash{$asmbl_id}->{'genomesize'}\t$ahref->{$asmbl_id}->{'gc'}\t$phageref->{$key}->{'5prime'}\t$phageref->{$key}->{'3prime'}\t$size\t$label\t$type\tN.D.\tN.D.\tN.D.\t$phageref->{$key}->{'gc'}\t$shortag\t";
##

  }
  $region2_header = "Best database match ($Best_hits{$best}->{'counts'}/$totalhits) to $dbref->{$best}->{'name'} \[$best\], taxonomy: $dbref->{$best}->{'taxon'}";
  $header = "$ahref->{$asmbl_id}->{'title'} [asmbl_id: $asmbl_id], $ahref->{$asmbl_id}->{'length'}, gc%: $ahref->{$asmbl_id}->{'gc'}%";
  $end5 = "";
  $end3 = "";
  $printfeatname = "";
  $annotation = "";
  write if ($DEBUG);
  write $FH;

  foreach $yek (sort {$a <=> $b} keys %{$phageref->{$key}->{'memberhash'}}) {
      print "$phage_end5 [ $yek ] $phage_end3\n" if ($DEBUG);
      $featname = $phageref->{$key}->{'memberhash'}{$yek}->{'featname'};
      if ((($yek >= $phage_end5) && ($yek <= $phage_end3)) || ($featname eq $phageref->{$key}->{'target'})) { # if gene is within phage boundaries +/- 1 gene or is the target
	  $end5 = $yek;
	  $end3 = $hitref->{$phage_hash{$key}->{'memberhash'}{$end5}->{'featname'}}->{'end3'};
	  if (($asmbly_status) && ($featname ne $phageref->{$key}->{'target'}) && ($phageref->{$key}->{'type'} eq "prophage"))  { # don't pull seqs for target
	      $seq = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $end5, $end3);
	      if (exists $hitref->{$featname}->{'genome_tag'}) { # new code 4/6/2010 so we have a clean feat_name-tag printed (ie giBLA-ABPI00000000_PFPR01)
		  $seq_header = "gi" . $featname . "-" . $hitref->{$featname}->{'genome_tag'} . "_" . $asmbl_id . "_" . $region_tag . " " . $hitref->{$featname}->{'com_name'}; # combodb format
	      }
	      else {
		  $seq_header = $featname . "_" . $asmbl_id . "_" . $region_tag . " " . $hitref->{$featname}->{'com_name'};
	      }
	      if (($phageref->{$key}->{'memberhash'}{$yek}->{'featname'} !~ /^att/) && ($featname !~ /RNA/)) { #skip over non-gene features and targets for integration
		  $phageref->{$key}->{'ORFcnt'}++; # count the number of ORFs per phage region
		  &printFasta(*SEQFILE, $seq_header, $seq);
		  ($pep) = &dna2pep($seq);
		  &printFasta(*PEPFILE, $seq_header, $pep);
	      }
	  }
	  if (exists $hitref->{$featname}->{'annotation'}) {
	      if ($featname =~ /^att/)  {  # don't print E-VALUE if an att site
		  $printfeatname = $phageref->{$key}->{'memberhash'}{$yek}->{'featname'};
		  $annotation = "$hitref->{$featname}->{'annotation'}";
		  write if ($DEBUG);
		  write $FH;
		  next;
	      }
	      else  {
		  $annotation = "$hitref->{$featname}->{'annotation'}, TAG = $hitref->{$featname}->{'phage'}, E-VALUE = $hitref->{$featname}->{'evalue'}";
	      }
	      $printfeatname = $featname;
	      write if ($DEBUG);
	      write $FH;
	      $end5 = "";
	      $annotation = "[ANNO] $hitref->{$featname}->{'com_name'}";
	      $printfeatname = "";
	      write if ($DEBUG);
	      write $FH;
	  }
	  else {
	      if ($featname =~ /tRNA/) {
		  my $com_feat = $asmbl_id . "_" . $featname;
		  $annotation = "[TRNA] $hitref->{$com_feat}->{'com_name'} $hitref->{$com_feat}->{'organism'}";
	      }
	      elsif ($featname =~ /tmRNA/) {
		  my $com_feat = $asmbl_id . "_" . $featname;
		  $annotation = "[TMRNA] $hitref->{$com_feat}->{'com_name'} $hitref->{$com_feat}->{'organism'}";
	      }
	      else {
		  $annotation = "[ANNO] $hitref->{$featname}->{'com_name'}";
	      }
	      $printfeatname = $featname;
	      write if ($DEBUG);
	      write $FH;
	  }
##
	  if (($phage_end5 <= $yek ) && ($print5 == 0) && (exists $phageref->{$key}->{'memberhash'}{$yek})) {
	      print TAB "$phageref->{$key}->{'memberhash'}{$yek}->{'featname'}\t";
	      $print5++;
	  }
	  elsif ($yek < $phage_end3) {
	      $hold_locus = $yek;
	  }
##
	  if ($featname eq $phageref->{$key}->{'target'}) { # if the gene is the target, note this!
	      $end5 = "";
	      $printfeatname = "";
	      $annotation = "[TARGET]";
	      write if ($DEBUG);
	      write $FH;
	  }
	  if (exists $hitref->{$featname}->{'hmm'}) {
	      foreach $HMM_name (sort {$a <=> $b} keys %{$hitref->{$featname}->{'hmm'}}) {
		  $hmm_hit_names = $hmm_hit_names . $HMM_name . "\t";
		  $end5 = "";
		  $printfeatname = "";
		  $annotation = "[HMM-$hitref->{$featname}->{'hmm'}{$HMM_name}->{'model'}] ($HMM_name): $hitref->{$featname}->{'hmm'}{$HMM_name}->{'hmm_com_name'}, SCORE = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'score'}, TRUSTED = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'trusted'}, NOISE = $hitref->{$featname}->{'hmm'}{$HMM_name}->{'noise'}";
		  write if ($DEBUG);
		  write $FH;
	      }
	  }
      }
  }
  print TAB "$phageref->{$key}->{'memberhash'}{$hold_locus}->{'featname'}\t$phageref->{$key}->{'int_HMM'}\t$phageref->{$key}->{'core_HMM'}\t$phageref->{$key}->{'above_noise_core_HMM'}\t$phageref->{$key}->{'lytic_HMM'}\t$phageref->{$key}->{'tail_HMM'}\t$phageref->{$key}->{'Mu_HMM'}\t$phage_direction\t$phageref->{$key}->{'att_distance'}\t$phageref->{$key}->{'ORFcnt'}\t$phageref->{$key}->{'serine_recomb'}\n";
##
};
if (($max > 1) || ($max == 0))  {
  print "There are $max putative large/medium/small phages in $asmbl_id!\n";
  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}
elsif ($max == 1)  {
  print "There is $max putative large/medium/small phage in $asmbl_id!\n";
  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

format OUTPUT_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format OUTPUT =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.
format OUTmedium_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.

format OUTmedium =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.

format OUTsmall_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format OUTsmall =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.

format STDOUT_TOP =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region_header
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$region2_header

END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: @<<
                                                                                                                                                                         $%
............................................................................................................................................................................
.
format STDOUT =
@<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<...
$end5,   $printfeatname,         $annotation
.
  foreach $key (sort {$a <=> $b} keys %{$phageref}) {        
      $label = "";
      $print_files = undef;
      $size = abs($phageref->{$key}->{'5prime'} - $phageref->{$key}->{'3prime'})+1;
      $size_att = abs($phageref->{$key}->{'5prime_att'} - $phageref->{$key}->{'3prime_att'})+1;
      print "print_regions: phage $key, size_att = $size_att\n" if ($DEBUG);
      if ($size_att > $size)  {  # if there is an att site and the size of the phage with att sites is bigger, use that to group small, med, large
	  $sizer = $size_att;
      }
      else  {
	  $sizer = $size;
      }
##
      if (exists $phageref->{$key}->{'direction'})  {
          $phage_direction = $phageref->{$key}->{'direction'};
      }
      else  {
          $phage_direction = "+"; # set default direction as -------> if no att found
      }
      if ($size_att > 1)  {
          if ($phage_direction eq "+") { # phage is ----->
              $phage_end5 = $phageref->{$key}->{'5prime_att'};
              $phage_end3 = $phageref->{$key}->{'3prime_att'};
              print "making + phage $key\n" if ($DEBUG);
          }
          else { # phage is <-----, so we want to switch phage_end5 to be the 5' most coord relative to the assembly
              $phage_end5 = $phageref->{$key}->{'3prime_att'};
              $phage_end3 = $phageref->{$key}->{'5prime_att'};
          }
          $phageref->{$key}->{'seq'} = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $phageref->{$key}->{'5prime_att'}, $phageref->{$key}->{'3prime_att'});
      }
      else  { # no att site found, use the coords in 5prime and 3prime
          if ($phage_direction eq "+") { # phage is ----->
              $phage_end5 = $phageref->{$key}->{'5prime'};
              $phage_end3 = $phageref->{$key}->{'3prime'};
	      print "Print_regions+:  >NO ATT< phage_end5 = ($phage_end5), phage_end3 = ($phage_end3)\n" if ($DEBUG);
          }
          else { # phage is <-----, so we want to switch phage_end5 to be the 5' most coord relative to the assembly
              $phage_end5 = $phageref->{$key}->{'3prime'};
              $phage_end3 = $phageref->{$key}->{'5prime'};
	      print "Print_regions-:  >NO ATT< phage_end5 = ($phage_end5), phage_end3 = ($phage_end3)\n" if ($DEBUG);
          }
          $phageref->{$key}->{'seq'} = &quickcut($ahref->{$asmbl_id}->{'sequence'}, $phageref->{$key}->{'5prime'}, $phageref->{$key}->{'3prime'});
      }

      if ($asmbly_status) {
        $phageref->{$key}->{'gc'} = &get_gc_content($phageref->{$key}->{'seq'});
      }
      else {
        $phageref->{$key}->{'gc'} = "NA";
      }
##

###
# check for type and whether has att or not
# 
     if (($size_att > 1) && (!defined ($phageref->{$key}->{'type'})))  {  # probable integrated element (plasmid, transposon, other)
           if ($phageref->{$key}->{'above_noise_core_HMM'} > 1) { #need 2 >= noise core HMM hits
              ## maybe make a more strict core set of portal, major capsid, terminase, capsid protease >= noise definition ##
              $phageref->{$key}->{'type'} = "prophage";
           }
           else  {
              $phageref->{$key}->{'type'} = "integrated element";
           }
      }
      elsif (($size_att > 1) && ($phageref->{$key}->{'type'} eq "prophage") && (!defined ($phageref->{$key}->{'class'})) && ($phageref->{$key}->{'core_HMM'} == 0) && ($phageref->{$key}->{'core_BLAST'} == 0))  {
           $phageref->{$key}->{'class'} = "satellite"; # if no core HMMs or BLAST hits are found, is a prophage with no class and an att site, call satellite (I know there are always exceptions)
      }
      elsif (($size_att == 1) && (!defined ($phageref->{$key}->{'type'}))) {
           $phageref->{$key}->{'type'} = "degenerate";
      }
###
      if ($key < 10) {
           $region_tag = "PFPR0$key";
      }
      else {
	   $region_tag = "PFPR$key";
      }
      $seq_tag = "PFPROPHAGE" . "_" . $bugtag . "_" . "$region_tag";
###

      if (($sizer < 10000) && ($strict == 0))  {
        if ($small_open == 0)  { # only open each file if needed and only once per invocation of the script so we can get multiple phages per file
          open (OUTsmall, ">$write_dir/$asmbl_id\_$hitsperwindow\_phpico\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phpico\_$window\_$step\_$evalue.out\n";
          $small_open = 1;
        }
	if ($sizer >= 5000) { $micro++; $totalmicrobp = $totalmicrobp + $sizer;}
        elsif ($sizer >= 2000) { $nano++; $totalnanobp = $totalnanobp + $sizer;}
        else { $pico++; $totalpicobp = $totalpicobp + $sizer;}
	$small++;
        $print_files = 1;
        $totalsmallbp = $totalsmallbp + $sizer;
        if (!defined ($phageref->{$key}->{'class'}))  {
          $label = "Small";
        }
        else {
          $label = "Small $phageref->{$key}->{'class'}";
        }
        &print_sub($key, \*OUTsmall, $label, $small, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
      }
      elsif (($sizer >= 10000) && ($sizer <= 18000))  {
        if (($strict == 0) || (($strict == 1) && ($phageref->{$key}->{'type'} eq "prophage")))  {
	  if ($med_open == 0)  {
            open (OUTmedium, ">$write_dir/$asmbl_id\_$hitsperwindow\_phmedio\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phmedio\_$window\_$step\_$evalue.out\n";
            $med_open = 1;
          }
	  $medium++;
          $print_files = 1;
          $totalmediumbp = $totalmediumbp + $sizer;
          if (!defined ($phageref->{$key}->{'class'}))  {
            $label = "Medium";
          }
          else {
            $label = "Medium $phageref->{$key}->{'class'}";
          }
          &print_sub($key, \*OUTmedium, $label, $medium, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
        }
        else {
	    print "Region $key is unlikely a prophage, deleting ...\n";
	    $phage_regions--;  # deprecate phage_region counter
	}
      }
      elsif ($sizer > 18000) { # Print the data for phage regions > 18 Kbp in size to a file, bug fix was reporting small regions as large
        if (($strict == 0) || (($strict == 1) && ($phageref->{$key}->{'type'} eq "prophage")))  {
	  if ($large_open == 0)  {
            open (OUTPUT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phregions\_$window\_$step\_$evalue.out") || die "can't open file $write_dir/$asmbl_id\_$hitsperwindow\_phregions\_$window\_$step\_$evalue.out\n";             
            $large_open = 1;
          }
          $print_files = 1;
	  $large++;
          $totalphagebp = $totalphagebp + $sizer;
          if (!defined ($phageref->{$key}->{'class'}))  {
            $label = "Large";
          }
          else {
            $label = "Large $phageref->{$key}->{'class'}";
          }
          &print_sub($key, \*OUTPUT, $label, $large, $phageref->{$key}->{'type'}, $phage_end5, $phage_end3, $phage_direction);
        }
	else { 
	    print "Region $key is unlikely a prophage, deleting ...\n";
	    $phage_regions--; # deprecate phage_region counter
	}
      }
      else {
	  print "Region $key is unlikely a prophage, deleting ...\n";
	  $phage_regions--;  # deprecate phage_region counter
      }
  if (($print_files) && ($phageref->{$key}->{'type'} eq "prophage"))  {
      $con_header = $ahref->{$asmbl_id}->{'title'};
      $con_header =~ s/\,.*$//; # remove everything past comma
      if (length($label) > 1) { $label = "$label "; }
      $con_header = $seq_tag . " " . $con_header . " " . $label . $phageref->{$key}->{'type'} . " region " . $key . " (" . $phage_end5 . "-" . $phage_end3 . " bp)";
      &printFasta(*CONFILE, $con_header, $phageref->{$key}->{'seq'});
  }
  }
     
$percent = ($totalphagebp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;
$percentsmall = ($totalsmallbp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;
$percentmedium = ($totalmediumbp/$asmbl_idhash{$asmbl_id}->{'genomesize'})*100;

print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

if ($phage_regions > "0")  {

    printf "There are $large regions > 18 Kb and summing $totalphagebp bp of sequence (%5.2f", $percent;
    print "% of the genome)\n";

    printf "There are $medium regions between 10 and 18 Kb and summing $totalmediumbp bp of sequence (%5.2f", $percentmedium;
    print "% of the genome)\n";

    printf "There are $small regions < 10 Kb and summing $totalsmallbp bp of sequence (%5.2f", $percentsmall;
    print "% of the genome)\n";

    if ($large > 0) {
	print OUTPUT "............................................................................................................................................................................\n";
	printf OUTPUT "There are $large regions > 18 Kb and summing $totalphagebp bp of sequence (%5.2f", $percent;
	print OUTPUT "% of the genome)\n";
	close (OUTPUT);
    }

    if ($medium > 0) {
	print OUTmedium "............................................................................................................................................................................\n";
	printf OUTmedium "There are $medium regions between 10 and 18 Kb and summing $totalmediumbp bp of sequence (%5.2f", $percentmedium;
	print OUTmedium "% of the genome)\n";
	close (OUTmedium);
    }

    if ($small > 0) {
	print OUTsmall "............................................................................................................................................................................\n";
	printf OUTsmall "There are $small regions < 10 Kb and summing $totalsmallbp bp of sequence (%5.2f", $percentsmall;
	print OUTsmall "% of the genome)\n";
	close (OUTsmall);
    }
}
}

sub write_output {  # write xgraph file and report file

  my ($write_dir,$asmbl_id,$hitsperwindow,$ahref,$figref) = @_;
  my $minimum_hits = $hitsperwindow-1;
  my $n = 0;
  my $max_peak = "";
  my $counts = 0;
  my $title = $ahref->{$asmbl_id}->{'title'};
  $title =~ s/$asmbl_id//;
  my $fig = "";
  my @feat_names = ();
  my $clean_counts = 0;

  open (OUTPUT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phgraph\_$window\_$step\_$evalue.out") || die "can't open file phgraph$window.out.\n";
  print OUTPUT "title = $title id $asmbl_id (>= $minimum_hits hits per window)\n";
  print OUTPUT "title_x = Position (bp)\n";
  print OUTPUT "title_y = # Hits per $window bp window/$step bp step\n";
  print OUTPUT "color = red\n";
  print OUTPUT "thickness = 1\n";
  open (REPORT, ">$write_dir/$asmbl_id\_$hitsperwindow\_phreport\_$window\_$step\_$evalue.out") || die "can't open file phreport$window.out.\n";
  foreach $fig (sort {$a <=> $b} keys %{$figref}) {
      $counts = $figref->{$fig}->{'counts'};
      if ($counts <= $minimum_hits)  {  # eliminates background hits of less than $hitsperwindow from graph
	  $clean_counts = 0;
      }
      else {
	  $clean_counts = $counts;
      }
      if (exists($figref->{$fig}->{'peak_number'})) {
	  $max_peak = $figref->{$fig}->{'peak_value'} + .1;
	  print OUTPUT "anno\t$fig\t$max_peak\t$figref->{$fig}->{'peak_number'}\n";
      }
      print OUTPUT "$fig\t$clean_counts\n";
      # strange syntax because a runtime exception occurs under use strict because you're dereferencing an undefined reference where autovivification won't occur.
      @feat_names = exists( $figref->{$fig}->{'featnames'} ) 
	  ? @{ $figref->{$fig}->{'featnames'} }
          : ();
      print REPORT "$fig\t$counts\t@feat_names\n";
  }
  close (OUTPUT);
  close (REPORT);
}

sub option_help {
   print <<_EOB_;
$prog - Find prophage regions in DNA sequence from btab (blastp) and hmm data
 Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved

 License:  This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Citation: Fouts, D. E. (2006) "Phage_Finder: automated identification and classification of prophage 
           regions in complete bacterial genome sequences." Nucleic Acids Res 34:5839-51. PMCID: PMC1635311.

 Usage: $prog <options>
 Example: $prog -t ncbi.out -i phage_finder_info.txt -r tRNAscan.out -n tmRNA_aragorn.out -A NC_000913.con -S
 Switch: -h for help\
 Option:
     -b: base directory path [default = PWD]
     -p: path to btab file (default = base directory)
     -t: name of WU-BLAST or NCBI (-m 8 option) btab input file [REQUIRED]
     -i: tab-delimitied flat file containing scaffold/contig/assembly_ID size_of_molecule feat_name end5 end3 com_name [REQUIRED]
     -m: htab file containing HMM data (REQUIRED for finding integrases and att sites)
     -F: search method (B or b for NCBI BLAST, M or m for MUMmer, F or f for FASTA33) (default = BLAST)
     -r: tRNAscan-SE output file [optional]
     -n: Aragon tmRNA-finding output file (-m option in aragon) [optional]
     -w: Scanning WINDOW size (default = 10000 nucleotides)
     -s: STEP size (default = 5000 nucleotides)
     -E: E-value (default = 0.000001)
     -H: Number of allowable hits per window to mark a region (default = 4)
     -a: User-defined asmbl_id to search (default picks asmbl_id with largest size)
     -A: File name of .1con
     -B: Path to .1con if not in base directory 
     -V: print version information
     -S: Strict mode:  print only regions that have core HMM hits or Mu-like and are > 10 Kbp (default = 0)
     -d: DEBUG MODE (default = 0)
 Format: 
 Output:
_EOB_
    exit;
}
################################## M A I N ####################################

print "$prog (October 2011) - identify prophage regions within bacterial genomes\n\n";
print "Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved\n";
print "This program comes with ABSOLUTELY NO WARRANTY; for details refer to LICENSE.txt.\n";
print "This is free software, and you are welcome to redistribute it\n";
print "under certain conditions; see LICENSE.txt for details.\n\n";
print "  Please cite:\n";
print "       Fouts, D. E. (2006) Phage_Finder: automated identification and classification\n"; 
print "       of prophage regions in complete bacterial genome sequences.\n";
print "       Nucleic Acids Res 34:5839-51. PMCID: PMC1635311.\n\n";

### Create log file in base directory ###
my $current_logfile = &create_log($logfile);
&write_log("0", $invocation);
### gather all information ###

print "Checking the format of $infofile ...\n" if ($DEBUG);
&write_log("1","Checking the format of $infofile");
$infofile = &check_infofile($basedir,$infofile,$asmbly_file,$defined_asmbl_id,$DEBUG);

print "Retrieving information from flat file $infofile . . .\n" if ($DEBUG);
&write_log("1","Retrieving information from flat file $infofile");
$num_contigs = &get_gene_info($infofile,$defined_asmbl_id,\%asmbl_idhash,\%hithash,\%rehash,$DEBUG);

print "There are $num_contigs contigs/assemblies in $infofile\n" if ($DEBUG);
&write_log("1","There are $num_contigs contigs/assemblies in $infofile");

if (defined($tRNA_data))  { # is there a tRNAscan-SE file?
  print "Looking for tRNA information ...\n" if ($DEBUG);
  &write_log("1","Looking for tRNA information");
  &get_tRNAs($tRNA_data,\%tRNA_hash,\%hithash,\%rehash,$DEBUG);
}
else {
  print "User did not specifiy a tRNAScan file or it contains no data, skipping any searches involving tRNA ...\n" if ($DEBUG);
  &write_log("1","User did not specifiy a tRNAScan file or it contains no data, skipping any searches involving tRNA");
}

if (defined($tmRNA_data))  { # is there an Aragorn file?
  print "Looking for tmRNA information ...\n" if ($DEBUG);
  &write_log("1","Looking for tmRNA information");
  &get_tmRNAs($tmRNA_data,\%tRNA_hash,\%hithash,\%rehash,$DEBUG);
}
else {
  print "User did not specifiy an Aragorn file or it contains no data, skipping any searches involving tmRNA ...\n" if ($DEBUG);
  &write_log("1","User did not specifiy an Aragorn file or it contains no data, skipping any searches involving tmRNA");
}

print "Getting list of OK phage common names ($comfile) . . .\n" if ($DEBUG);
&write_log("1","Getting list of OK phage common names ($comfile)");
&get_ok_comnames($comfile,\%ok_comnames,$DEBUG);

print "Getting tags, names and taxon info from BLAST DB file ($DBfile) . . .\n" if ($DEBUG);
&write_log("1","Getting tags, names and taxon info from BLAST DB file ($DBfile)");
&get_DB_info($DBfile,\%DB_info,$DEBUG);

print "Obtaining hits to phage database from $btabfile . . .\n" if ($DEBUG);
&write_log("1","Obtaining hits to phage database from $btabfile");
&select_featnames_from_btab($btabfile,$evalue,\%hithash,\%searchhash,\%exclude_hash,\%lists_hash,$DEBUG);

if (defined($hmmfile))  { # is there an hmm file?
    open (HMMFILE, "<$hmmfile") || &write_log("4","can't open file HMM datafile $hmmfile: $!\n");
       while (<HMMFILE>)  {
           if (/^# HMMER 3/) {
              $HMMversion = "3";
              last;
           }
           else {
              $HMMversion = "2";
              last;
           }
       }
    close (HMMFILE);
    push (@hmm_data, $hmmfile);
}
elsif (-s "$basedir/combined.hmm3") {
    $HMMversion = "3";
    push (@hmm_data, "$basedir/combined.hmm3");
}
else {
    $HMMversion = "2";
    push (@hmm_data, "$basedir/combined.hmm_GLOCAL");
    push (@hmm_data, "$basedir/combined.hmm_FRAG");
}

#### loop through each filename in the @filenames array in the load_hash_from_file subroutine

if ($HMMversion == "2") {
   @filenames = ("phage_exclude_v2.lst",
                   "core_hmm.lst",
                   "lysin_holin_hmm.lst",
                   "tails_hmm.lst",
                   "large_term_hmm.lst",
                   "small_term_hmm.lst",
                   "portal_hmm.lst",
                   "protease_hmm.lst",
                   "capsid_hmm.lst",
                   "Small_term_v2.lst",
                   "Large_term_v2.lst",
                   "Portal_v2.lst",
                   "Protease_v2.lst",
                   "Scaffold_v2.lst",
		   "Major_capsid_v2.lst");
}
elsif ($HMMversion == "3") {
   @filenames = ("phage_exclude_v2.lst",
                   "core_hmm3.lst",
                   "lysin_holin_hmm3.lst",
                   "tails_hmm3.lst",
                   "large_term_hmm3.lst",
                   "small_term_hmm3.lst",
                   "portal_hmm3.lst",
                   "protease_hmm3.lst",
                   "capsid_hmm3.lst",
                   "Small_term_v2.lst",
                   "Large_term_v2.lst",
                   "Portal_v2.lst",
                   "Protease_v2.lst",
                   "Scaffold_v2.lst",
                   "Major_capsid_v2.lst");
}

&load_hash_from_file($phome,\@filenames,\%exclude_hash,\%lists_hash,\%HMMs_hash,$DEBUG);

####

if (! &find_hmms($home,$phome,$HMMversion,\@hmm_data,$infofile,\%hithash,\%searchhash,\%HMMs_hash,\%serine_HMM_hash,\%tyrosine_HMM_hash,$DEBUG)) {
    print "User did not specify an hmm file used in finding integrases, skipping this analysis . . .\n" if ($DEBUG);
    &write_log("1","User did not specify an hmm file used in finding integrases, skipping this analysis");
  }

### ---------------------------- ###

# all these datafiles were pulled from the contig directories into the main genome directory to facilitate access to this data when there were many contigs (ie in draft genomes)
# their systematic names will make it easier for a bash script to concatenate all of them from a listfile

my $apisfile = "$basedir/PFPR.apis";
my $tabfile = "$basedir/PFPR_tab.txt"; # added 08/22/2008 by Derrick Fouts
my $smalltermpep = "$basedir/PFPR_small_terminase.pep";
my $largetermpep = "$basedir/PFPR_large_terminase.pep";
my $portalpep = "$basedir/PFPR_portal.pep";
my $proteasepep = "$basedir/PFPR_protease.pep";
my $scaffpep = "$basedir/PFPR_scaffold.pep";
my $capsidpep = "$basedir/PFPR_capsid.pep";
my $serinerecombpep = "$basedir/PFPR_serine_recombinase.pep"; # only the one we think is used for integration is pulled out
my $tyrosinerecombpep = "$basedir/PFPR_tyrosine_recombinase.pep"; # only the one we think is used for integration is pulled out

open (TAB, ">$tabfile") || die "can not open file $tabfile for writing\n";
print TAB "#asmbl_id\tgenome_size\tgenome_gc\tbegin_region\tend_region\tsize_region\tlabel\ttype\t5prime_att\t3prime_att\ttarget\tregion_gc\tbest_db_match\tbegin_gene\tend_gene\t#integrase_HMMs\t#core_HMMs\t#>noise_HMMs\t#lytic_HMMs\t#tail_HMMs\t#Mu_HMMs\tregion_orientation\tdistance_int_to_att\t#genes\t#serine_recombinases\n"; # print header line for PFPR.tab file

open (CONFILE, ">PFPR.con") || die "can't open file PFPR.con\n";
open (SEQFILE, ">PFPR.seq") || die "can't open file PFPR.seq\n";
open (PEPFILE, ">PFPR.pep") || die "can't open file PFPR.pep\n";

### now, loop through assemblies, looking for phages ###
foreach $asmbl_id (sort {$a <=> $b} keys %asmbl_idhash)  {
  $num_contigs--;
  $write_dir = &create_dir($basedir,"$asmbl_id");
  ### main loop to count hits within window, incrementing by step size ### 
  %fighash = (); # clear out %fighash hash for each asmbl_id/contig
  %phage_hash = (); # clear our %phage_hash for each asmbl_id/contig
  @genomearray = (); # clear our @genomearray for each asmbl_id/contig
  for (my $n = 0; $n <= $asmbl_idhash{$asmbl_id}->{'genomesize'}; $n = $n + $step)  {  
    $fighash{$n}->{'counts'} = 0;
    foreach my $keys (sort {$a <=> $b} keys %{$searchhash{$asmbl_id}}) {
       if (($keys > $n) && ($keys < ($n+$window)))  {
         $fighash{$n}->{'counts'}++;
         $fighash{$n}->{'featnames'}[($fighash{$n}->{'counts'}-1)] = $searchhash{$asmbl_id}{$keys}->{'featname'};  # for each step, store the orf names with hits 
       }
       elsif ($keys >= ($n+$window)) { last; } # end foreach loop when keys (end5 is greater than or equal to window size
       if ($keys < ($n+$step)) { delete $searchhash{$asmbl_id}{$keys}; }  # Remove keys that have already been counted and are less than step size so not resorted and counted again  
    }
  }
  $hitsperwindow = $default_hpw;  # reset the user-defined or default hitsperwindow 
  if (($asmbl_idhash{$asmbl_id}->{'genomesize'} < 10000) && ($asmbl_idhash{$asmbl_id}->{'genomesize'} >= 7500))  {$hitsperwindow = 3;}

  $phage_regions = &find_regions($asmbl_id,$hitsperwindow,\%fighash,\%hithash,\%phage_hash,$DEBUG);

  if (exists ($phage_hash{1}->{'seed_end5'})) { # added 11/14/02 to not to attempt extension of phage regions when there are none! 
    &populate_genomearray(\%rehash,\@genomearray,\%hithash,$asmbl_id);
    &find_5prime_end($phage_regions,$asmbl_id,\@genomearray,\%phage_hash,\%hithash,\%tRNA_hash,\%rehash,\%ok_comnames,\%fighash,$DEBUG);  # find 5' boundary of phage regions
    &find_3prime_end($phage_regions,$asmbl_id,\@genomearray,\%phage_hash,\%hithash,\%tRNA_hash,\%rehash,\%ok_comnames,\%fighash,$DEBUG);  # find 3' boundary of phage regions
    &write_log("1","Writing output files");

    if (!defined($asmbly_file)) {
      print "WARNING: User did not specify a .1con file used in finding att sites, skipping this analysis . . .\n";
      &write_log("1","WARNING: User did not specify a .1con file used in finding att sites, skipping this analysis");
      $asmbly_status = 0;
    }
    else { $asmbly_status = 1 };

    $asmbly_status = &get_assemblies($asmbly_file,$infofile,\%assembly_hash,\%asmbl_idhash,$asmbl_id,$DEBUG) if ($asmbly_status == 1); # if we have a defined asmbl_id, then get the sequence...

    if ($asmbly_status == 1)  {
      print "..................................................................................................\n";
      print "Working on $assembly_hash{$asmbl_id}->{'title'} assembly|contig|scaffold id $asmbl_id which is $assembly_hash{$asmbl_id}->{'length'} bp in size, gc%: $assembly_hash{$asmbl_id}->{'gc'}%\n";
      &write_log("1","Working on $assembly_hash{$asmbl_id}->{'title'} assembly|contig|scaffold id $asmbl_id which is $assembly_hash{$asmbl_id}->{'length'} bp in size, gc%: $assembly_hash{$asmbl_id}->{'gc'}%");
    }
    &determine_region_type(\%phage_hash,\%hithash,\%HMMs_hash,\%serine_HMM_hash,\%tyrosine_HMM_hash,$DEBUG);
    if (((defined($hmmfile)) || (defined($tRNA_data))) && ($asmbly_status == 1))  { # must have either hmm data or tRNA data AND an assembly to work with to look for att sites
      print "===> ASMBLY_STATUS: $asmbly_status\n" if ($DEBUG);
      print "Looking for putative phage attachment (att) sites using $search_method . . .\n" if ($DEBUG);
      &write_log("1","Looking for putative phage attachment (att) sites using $search_method");
      $prefix = "$write_dir/$asmbl_id\_$hitsperwindow\_phtest\_$window\_$step\_$evalue";
      &find_att_sites($prefix,$asmbl_id,$search_method,\%assembly_hash,\%phage_hash,\%hithash,\%rehash,\%tRNA_hash,\@genomearray,\%fighash,$window,$step,$hitsperwindow,\%ok_comnames,$DEBUG);
    }
    &print_regions($write_dir,$asmbly_status,$asmbl_id,$hitsperwindow,$window,$step,$evalue,$phage_regions,\%phage_hash,\@genomearray,\%assembly_hash,\%hithash,\%rehash,\%DB_info,$DEBUG);    # print output of analysis
    &write_output($write_dir,$asmbl_id,$hitsperwindow,\%assembly_hash,\%fighash);     # print file for xgraph input
    &write_log("2", $num_contigs, $asmbl_id);
  }
##############
  if ($phage_regions == "0")  {  # added 11/14/02 to give default output of NO phages if there are no predicted prophages (ie ntcj01 has none)
    print "Sorry, no phages in $asmbl_id :( . . .\n";
    unlink glob "$write_dir/*";  # remove the directory if no prophages identified (to keep it clean for metagenomic studies)
    rmdir $write_dir;
    &write_log("3", $num_contigs, $asmbl_id);
  }
delete $assembly_hash{$asmbl_id}; # remove old assembly information (sequence, title, ...)
##############
}
close (TAB);
close (CONFILE);
close (SEQFILE);
close (PEPFILE);

if (-z "PFPR.con" == 1)  { # remove empty files
    print "removing PFPR.con ...\n" if ($DEBUG);
    unlink "PFPR.con";
}
if (-z "PFPR.seq" == 1)  {
    print "removing PFPR.seq ...\n" if ($DEBUG);
    unlink "PFPR.seq";
} 
if (-z "PFPR.pep" == 1)  {
    print "removing PFPR.pep ...\n" if ($DEBUG);
    unlink "PFPR.pep";
}
exit(0);

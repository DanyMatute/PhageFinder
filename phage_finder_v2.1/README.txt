Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

CREDITS
-------
Math::Round was written by Geoffrey Rommel and is available at:
(http://search.cpan.org/~grommel/Math-Round-0.06/Round.pm).

PATCHES
-------

Suggested updates should be directed to Derrick Fouts (dfouts@jcvi.org) for consideration.
04/09/2012:  Fixed a bug in the Perl module Phage_subs_v2.pm that prevented the proper creation and usage of the phage_finder_info.txt file when providing a Genbank .ptt file.
          :  Corrected HMM3_searches.pl (hmm3search changed to hmmsearch).
          :  phage_finder_v2.1.sh now checks for correct hmmsearch version.  Errors out if not "HMMER 3"

INTRODUCTION
------------

Phage_Finder is a heuristic computer program written in PERL that uses BLASTP data, HMM results, and tRNA/tmRNA information to find prophage regions in complete bacterial genome sequences.  For more information please visit the Phage_Finder website at http://phage-finder.sourceforge.net/

Phage_Finder was written by: 

Derrick E. Fouts, Ph.D.
Genomic Medicine
The J. Craig Venter Institute (JCVI)
9704 Medical Center Drive
Rockville, MD  20850
(301) 795-7874
dfouts@jcvi.org

SYSTEM REQUIREMENTS
-------------------

The programs should run on all Unix platforms.  It has been tested on Linux (CentOS/Redhat, Suse and Debian) and Mac OS X 10.3-10.6 operating systems.  It has not been tested on all Unix platforms.

SOFTWARE REQUIREMENTS/DEPENDENCIES
----------------------------------

Phage_Finder_v2.1.pl requires the following programs or packages for full functionality:

PERL version 5.8.5 or later (http://www.perl.org)
NCBI BLASTALL version 2.2.10 or later [Linux] (ftp://ftp.ncbi.nih.gov/blast/executables/release/) or [Mac OS X] Apple/Genentech optimized NCBI BLAST (http://developer.apple.com/opensource/tools/blast.html) 
WUBLAST 2.0 (now AB-BLAST at http://blast.advbiocomp.com/)
HMMSEARCH (http://hmmer.janelia.org/) (this is the new HMMER3)
tRNAscan-SE (http://selab.janelia.org/tRNAscan-SE/)
Aragorn version 1.2 (http://130.235.46.10/ARAGORN/)
FASTA33, MUMMER if you want to use these to find attachment sites
Math::Round PERL module (http://search.cpan.org/~grommel/Math-Round-0.06/Round.pm) [version 0.05 INCLUDED in this distribution]
PHAGE::Phage_subs_v2 PERL module [INCLUDED in this distribution]
Getopt::Std [should be installed with PERL]
XGRAPH (http://www.xgraph.org/)

INCLUDED IN DISTRIBUTION
------------------------

-PERL SCRIPTS AND MODULES- 

~/phage_finder_v2.1/bin/Phage_Finder_v2.1.pl:  The main PERL script for finding prophage regions

~/phage_finder_v2.1/lib/PHAGE/Phage_subs_v2.pm:  PERL module that contains several reusable subroutines need by Phage_Finder_v2.1.pl (NOT object oriented)

~/phage_finder_v2.1/lib/Math/Round.pm:  The Math::Round PERL module by Geoffrey Rommel.  Phage_Finder uses the nhimult and nlowmult to round numbers up or down by a defined multiple.  Very cool module Geoffrey :)

-BASH SHELL SCRIPTS-

~/phage_finder_v2.1/bin/phage_finder_v2.1.sh:  A BASH shell script used to run the entire Phage_Finder pipeline (including BLAST, HMM-finding, tRNA/tmRNA-finding and the Phage_Finder_v2.1.pl script itself)

~/phage_finder_v2.1/bin/HMM3_searchs.sh:  A BASH shell script to run all of the phage HMMs using HMMER3, reporting the progress as % completed and concatenating the results into a combined.hmm3 file

-BLAST DATABASE-

~/phage_finder_v2.1/DB/phage_10_02_07_release.db:  NCBI BLAST-formatted (.phr, .pin, and .psq) and WUBLAST-formatted (.ahd, .atb, and .bsq) files.

-Helper ascii files-

~/phage_finder_v2.1/hmm3.lst:  File containing a list of all 322 HMM3 phage HMMs.  14 old HMM2 were merged with existing and new HMM2 models

~/phage_finder_v2.1/HMM3_master.lst:  File containing the trusted and noise cutoff values and the name of each phage HMMER3 HMM

~/phage_finder_v2.1/phage_com_names_combo.txt:  List of acceptable phage annotation

~/phage_finder_v2.1/phage_exclude_v2.list:  List of accessions to exclude from analysis when using the 10_02_07 curated BLAST database

~/phage_finder_v2.1/PHAGE_core_HMM3.lst:  List of "core" phage HMMMER3 HMMs

~/phage_finder_v2.1/lysin_holin_hmm3.lst:  List of lysis or holin phage HMMER3 HMMs

~/phage_finder_v2.1/tails_hmm3.lst:  List of phage tails HMMER3 HMMs

~/phage_finder_v2.1/terminase_hmm3.lst:  List of phage large terminase HMMER3 HMMs

~/phage_finder_v2.1/portal_hmm3.lst:  List of phage portal HMMER3 HMMs

~/phage_finder_v2.1/large_term_hmm3.lst:  List of phage large terminus HMMER3 HMMs

- Note - All lists of manually-curated Genbank accessions are provisional with no guarantees of correctness.  This is a work in progress.

~/phage_finder_v2.1/Large_term_v2.lst:  List of manually-curated phage Large terminase accessions for 10_02_07 db

~/phage_finder_v2.1/Major_capsid_v2.lst:  List of manually-curated phage Major capsid accessions for 10_02_07 db

~/phage_finder_v2.1/Portal_v2.lst:  List of manually-curated phage portal accessions for 10_02_07 db

~/phage_finder_v2.1/Protease_v2.lst:  List of manually-curated phage capsid protease accessions for 10_02_07 db

~/phage_finger_v2.0/Scaffold_v2.lst:  List of manually-curated phage scaffold accessions for 10_02_07 db

~/phage_finder_v2.1/Small_term_v2.lst:  List of manually-curated phage Small terminus accessions for 10_02_07 db

-Hidden Markov Models (HMMs)-

~/phage_finder_v2.1/PHAGE_HMM3s_dir/:  Directory containing 322 phage HMMER3 HMM models

~/phage_finder_v2.1/examples_dir/:  Directory containing the legacy 42-genome test dataset

INSTALLATION
------------

First, place the distribution tarball to your home directory (~/)

Second, uncompress the distribution tarball by typing:

% tar -xvzf phage_finder_v2.tar.gz

REQUIRED INPUT FILES
--------------------

1) WU-BLAST or NCBI (-m 8 option) btab input file
2) phage_finder_info.txt file (a tab-delimitied file containing scaffold/contig/assembly_ID size_of_molecule feat_name end5 end3 com_name)

INVOCATION
----------

 Usage: Phage_Finder_v2.1.pl <options>
 Example: Phage_Finder_v2.1.pl -t ncbi.out -i phage_finder_info.txt -r tRNAscan.out -n tmRNA_aragorn.out -A NC_000913.con -S

 Switch: -h for help
 Options:
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
 Output:  All stored within a subdirectory of the current working directory ($PWD) named by the genome contig or accession id (ie NC_000913)
          1) phage_phinder_<id>.log:  a log file recording Phage_Finder progress
          2) phgraph file:  an XGRAPH plot of the phage regions
          3) phreport file:  a tab-delimited report file that shows (coordinate incremented by the step size, # hits per window, and the feat_name or locus name of the hits) 
          4) phpico, phmedio, phregions:  tab-delimited files containing the 5’ end of each gene, tRNA or att site within each region, the name of the feature, and the annotation/database match/HMM match as well as the G+C% content of each region, a best guess for the type of region, and the coordinates of each region with or without att site adjustments. There are three different names for this file, depending on the size of the regions (1-10000 bp [phpico], 10001-18000 bp [phmedio] and >18001 bp [phreigons])
          5) PFPR_tab.txt file:  a tab-delimited file containing (contig_id, size of the genome, G+C% content of the genome, begin coord. of the phage region, end coord. of the phage region, size of region in bp, label (small, medium, large), region type (prophage, integrated element, degenerate), sequence of attR, sequence of attL, name of integration target, G+C% of region, best db match (NEW), begin feat_name or locus name, end feat_name or locus name, # integrase HMM hits, # core_HMM hits, # above noise core_HMM hits, # lytic gene HMM hits, # tail HMM hits, # Mu HMM hits, orientation of the prophage based on orientation of the target or the position of the integrase, the distance from att site to integrase, the number of genes in the region, and number of serene recombinases (NEW)
          6) PFPR.con file:  a file in FASTA format containing the DNA sequence of all phage regions
          7) PFPR.seq file:  a file in FASTA format containing the DNA sequence of each gene with each phage region
          8) PFPR.pep file:  a file in FASTA format containing the protein sequence of each gene within the phage region

PHAGE_FINDER PIPELINE
---------------------

To run the complete Phage_Finder pipeline automaticly, you need to have the following:

1) each genome within a directory that is named by accession or contig_id (ie NC_000913) (see examples_dir)
2) protein sequences for each genome in their repective directories with extension (.pep or .faa)
3) phage_finder_info.txt or GenBank .ptt file in their repective directories
4) complete genome sequence (.con or .fna file) in their repective directories

% Phage_Finder_v2.1.sh <file with list of accessions>

EXAMPLE DATA
------------

Sample data can be found in ~/phage_finder_v2.1/examples_dir

The file phage_phinder_postprocess.out contains data parsed using PERL from each test genome's PFPR_tab.txt file generated by Phage_Finder.

Each column (number by PERL convention) in this file represents:

0) ORGANISM = The name of the organism tested
1) ACCESSION = The GenBank accession number of the test genome
2) SIZE = The size in base pairs of the test genome
3) GC% = The G+C% of the test genome
4) # = The number of the prophage region
5) pgc% = The G+C% of the prophage region
6) end5 = The 5' end of the prophage region
7) end3 = The 3' end of the prophage region
8) loci = The gene span of the prophage region (locus 5'-locus 3')
9) ori = The orientation of the prophage region
10) #phgbp = The size (in base pairs) of the prophage region
11) Ty = Region type (PRO=prophage, BAC=bacteriocin, DEG=degenerate, IEL=integrated element)
12) class = Region class (LRG=large, MED=medium, SML=small, SAT=satellite, RET=retron element, Mu=Mu-like, P2=P2-like, P4=P4-like)
13) att = Predicted attachment site found (Y/N)
14) target = Predicted target of insertion
15) In = The number of integrase HMM hits
16) Co = The number of core phage HMM hits above the trusted cutoff
17) NC = The number of core phage HMM hits above the noise cutoff
18) Ly = The number of phage lysis HMM hits
19) Ta = The number of phage tail HMM hits
20) Mu = The number of Mu-like phage HMM hits
21) dX = The distance from integrase to attachment site in base pairs
22) OC = The total number of open reading frames (ORFs) within region

You can also see the legacy results from searching 302 complete bacterial genomes in the file ~/phage_finder_v2.1/examples_dir/302_bac_genomes.txt


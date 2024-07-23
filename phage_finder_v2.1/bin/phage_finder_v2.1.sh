#!/bin/sh

# Bash script to run the Phage_Finder pipeline with HMMer3 models

# Usage: phage_finder_v2.0.sh <prefix of .pep/.ffa, .ptt and .con/.fna file>

# NOTE: a phage_finder_info.txt file will be searched before a .ptt file
# .pep is the multifasta protein sequence file
# .ptt is a GenBank .ptt file that has the coordinates and ORF names with annotation
# .con is a file that contains the complete nucleotide sequence of the genome being searched

home=`echo $HOME`
ver="2.1"
phome="$home/phage_finder_v$ver"
base=`pwd`
prefix=$1

  if [ -s $base/$prefix.pep ] # check if .pep file is present
  then
      pepfile="$prefix.pep"
  elif [ -s $base/$prefix.faa ]
  then
      pepfile="$prefix.faa"
  else
     echo "Could not file $prefix.pep or $prefix.faa.  Please check to make sure the file is present and contains data"
     exit 1
  fi 
    if [ -s $base/phage_finder_info.txt ] # check for phage_finder info file and if it has contents
    then
        infofile="phage_finder_info.txt"
    elif [ -s $base/$prefix.ptt ]
    then
          infofile="$prefix.ptt"
    else
      echo "Could not find a phage_finder_info.txt file or $prefix.ptt file.  Please make sure one of these files is present and contains data."
      exit 1
    fi
    if [ ! -s $base/combined.hmm3 ] # if GLOCAL HMM results not present, search
    then
        ## conduct GLOCAL HMMER3 searches
	HMMversion=`hmmsearch -h | head -2 | perl -ne 'chomp; if (/^#\s(HMMER\s\d)\./) {print "$1\n";'}`
	if [ "$HMMversion" != "HMMER 3" ]
	then
          echo "ERROR: This version of Phage_Finder only workes on HMM data generated using HMMER version 3.  Please check to make sure HMMER3 is installed and available using hmmsearch.  You can check the version by typing <hmmsearch -h | head -1>"
     	  exit 1
        fi 
        echo "  HMMER3 searches ..."
        $phome/bin/HMM3_searches.sh $base/$pepfile $ver
    fi
    if [ ! -e $base/ncbi.out ] # if BLAST results not present, search
    then
        ## do NCBI BLASTP searches
        echo "  BLASTing $pepfile against the Phage DB ..."
        blastall -p blastp -d $phome/DB/phage_10_02_07_release.db -m 8 -e 0.001 -i $pepfile -o ncbi.out -v 4 -b 4 -a 2 -F F
    fi
    if [ -s $base/$prefix.con ]
    then
        contigfile="$prefix.con"
    elif [ -s $base/$prefix.fna ]
    then
        contigfile="$prefix.fna"
    else
        echo "Could not find a phage_finder_info.txt file or $prefix.ptt file.  Please make sure one of these files is present and contains data.  In the meantime, I will go ahead and run phage_finder.pl without this information, but beware... NO att sites will be found!"
        contigfile=""
    fi
    if [ ! -e $base/tRNAscan.out ] && [ $base/$contigfile ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tRNAs
        echo "  find tRNA sequences ..."
        tRNAscan-SE -B -o tRNAscan.out $base/$contigfile > /dev/null
    fi

    if [ ! -e $base/tmRNA_aragorn.out ] && [ -e $base/$contigfile ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tmRNAs
        echo "  find tmRNA sequences ..."
        aragorn -m -o tmRNA_aragorn.out $base/$contigfile
    fi

    ## find the phage
    echo "  searching for Prophage regions ..."
    $phome/bin/Phage_Finder_v$ver.pl -t ncbi.out -i $infofile -r tRNAscan.out -n tmRNA_aragorn.out -A $contigfile -S
    rm error.log formatdb.log

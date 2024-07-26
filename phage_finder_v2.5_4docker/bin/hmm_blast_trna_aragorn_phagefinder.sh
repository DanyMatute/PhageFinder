#!/bin/sh
############################################################
# Help                                                     #
############################################################
Help()
{
    # Display Help
    echo "############################################################"
    echo "# HELP                                                     #"
    echo "############################################################"
    echo
    echo "DESCRIPTION: This script will accept sample.faa/pep, sample.fna/con, phage_finder_info.txt/ptt, tmRNA_aragorn.out, and and optional a tRNAscan.out file. The files will serve as inputs to conduct a HMMER3 search, NCBI BLASTP search, Aragron tmRNA search and an optional tRNAscan-SE tmRNA search. The output of the above commands will feed into Phage_Finder_v2.5.pl"
    echo
    echo
    echo "USAGE: hmm_blast_trna_aragorn_phagefinder.sh [-e] <sample.faa or sample.pep> <sample.fna or sample.con> <phage_finder_info.txt or sample.ptt> <tmRNA_aragorn.out>"
    echo 
    echo "OPTIONS:"
    echo "  -t, --trnascan tRNAscan.out   Use own tRNAscan.out, otherwise default will be used. "
    echo
    echo "############################################################"
}
############################################################
############################################################
# Main program                                             #
############################################################
############################################################
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":h" option; do
    case $option in
        h) # display Help
            Help
            exit;;
    esac
done

echo "Hello world!"
        
# # hmm searches
# echo "Running hmmsearch";
# for i in `cat /opt/PhageFinder/phage_finder_v2.5_4docker/hmm3.lst`; do
#     hmmsearch /opt/PhageFinder/phage_finder_v2.5_4docker/PHAGE_HMM3s_dir/$i.HMM {wildcards.genome}.faa;
# done > combined.hmm3;
        
# # blast
# echo "Running blastall";
# blastall -p blastp -d {pfHome}/DB/phage_10_02_07_release.db -m 8 -e 0.001 -i {wildcards.genome}.faa \
#     -o ncbi.out -v 4 -b 4 -a 2 -F F;

# # tRNA scan
# echo "Running tRNAscan-SE";
# if [[ -e tRNAscan.out ]]; then rm tRNAscan.out; fi
# tRNAscan-SE -B -Q -o tRNAscan.out {wildcards.genome}.fna;

# # aragorn
# echo "Running aragorn";
# aragorn -m -o tmRNA_aragorn.out {wildcards.genome}.fna;

# # phage_finder
# echo "Running phage_finder";
# {pfRun} -t ncbi.out -i phage_finder_info.txt -r tRNAscan.out -n tmRNA_aragorn.out -A {wildcards.genome}.fna -S;
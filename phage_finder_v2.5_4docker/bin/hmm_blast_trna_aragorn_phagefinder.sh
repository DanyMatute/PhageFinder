#!/bin/bash
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
    echo "USAGE: hmm_blast_trna_aragorn_phagefinder.sh [-t] <sample.faa or sample.pep> <sample.fna or sample.con> <phage_finder_info.txt or sample.ptt> <outputfile_location>"
    echo 
    echo "OPTIONS:"
    echo "  -t, --trnascan <tRNAscan.out>   Perform tRNAscan.out, otherwise default will be used. (optional) "
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
task_hmmsearch(){
    echo "  >>  Running hmmsearch"
    while IFS= read -r i; do
        echo "${i}, ${1}"
        hmmsearch /opt/PhageFinder/phage_finder_v2.5_4docker/PHAGE_HMM3s_dir/"$i".HMM "$1"
    done < /opt/PhageFinder/phage_finder_v2.5_4docker/hmm3.lst > combined.hmm3
    echo "  >>  hmmsearch DONE"
}

task_blast(){
    echo "  >>  Running BLASTp"
    blastp -db /opt/PhageFinder/phage_finder_v2.5_4docker/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query "$1" -out ncbi.out -max_target_seqs 5 -num_threads 8 
    echo "  >>  BLASTp DONE"
}

task_trnascan(){
    echo "  >>  Running tRNA-Scan"
    blastp -db /opt/PhageFinder/phage_finder_v2.5_4docker/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query "$1" -out ncbi.out -max_target_seqs 5 -num_threads 8 
    echo "  >>  tRNA-Scan DONE"
}

main() {
    # Positional Arguments
    local faa_file=
    local fna_file=
    local phagefinder_info_file=
    local outputfile_location=

    # Flag arguments
    local transcan_file="false"
    local position=0


    # Count all the args passed, and loop while we have all the args.
    while [[ "${#}" -gt 0 ]]; do   
        case "${1}" in 
            -h|--help|help)
                Help
                exit 0
                ;;
            -t|--trnascan)
                transcan_file="true"
                shift
                ;;
            
            *)
                case "${position}" in
                    0)
                        faa_file="${1}"
                        position=1
                        shift
                        ;;
                    1)
                        fna_file="${1}"
                        position=2
                        shift
                        ;;
                    2)
                        phagefinder_info_file="${1}"
                        position=3
                        shift
                        ;;
                    3)
                        outputfile_location="${1}"
                        position=4
                        shift
                        ;;
                    4) 
                        printf "    !! Unknown argument: %s\n\n" "${1}" >&2
                        Help >&2
                        exit 1
                        ;;
                esac
                ;;
        esac
    done

    # Validation

    [[ -z "${outputfile_location}" ]] && printf "    !!  Input Missing\n\n" >&2 && Help >&2 && exit 1

    # Echo inputs
    echo  "faa_file: ${faa_file}"
    echo  "fna_file: ${fna_file}"
    echo  "phagefinder_info_file: ${phagefinder_info_file}"
    echo  "outputfile_location: ${outputfile_location}"
    [[ "${transcan_file}" == "true" ]] && echo "transcan_file: ${transcan_file}"

    task_hmmsearch "${faa_file}"
    task_blast "${faa_file}"
    [[ "${transcan_file}" == "true" ]] && task_trnascan "${faa_file}"
    
    return 0

    

}

main "${@:-}"





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
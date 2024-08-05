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
    echo "DESCRIPTION: This script will accept sample.faa/pep, sample.fna/con, phage_finder_info.txt/ptt files and an output location. The files will serve as inputs to conduct a HMMER3 search, NCBI BLASTP search, Aragron tmRNA search and an optional tRNAscan-SE tmRNA search. The output of the above commands will feed into Phage_Finder_v2.5.pl"
    echo
    echo
    echo "USAGE: hmm_blast_trna_aragorn_phagefinder.sh [-t] <sample.faa or sample.pep>  <sample.fna or sample.con>  <phage_finder_info.txt or sample.ptt>  <outputfile_location>"

    echo 
    echo "OPTIONS:"
    echo "  -t, --trnascan      Perform tRNAscan-SE to produce tRNAscan.out and include in Phage_finder."
    echo
    echo "############################################################"
}

############################################################
############################################################
# Tasks                                           #
############################################################
############################################################

task_hmmsearch(){
    echo  ""
    echo "  RUNNING hmmsearch"
    while IFS= read -r i; do
        echo "${i}, ${1}"
        hmmsearch /opt/PhageFinder/phage_finder_v2.5_4docker/PHAGE_HMM3s_dir/"$i".HMM "$1"
    done < /opt/PhageFinder/phage_finder_v2.5_4docker/hmm3.lst > "${2}/combined.hmm3"
    echo "  > hmmsearch DONE"
}

task_blast(){
    echo  ""
    echo "  Running BLASTp"
    blastp -db /opt/PhageFinder/phage_finder_v2.5_4docker/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query "$1" -out "${2}/ncbi.out" -max_target_seqs 5 -num_threads 8 
    echo "  > BLASTp DONE"
}

task_trnascan(){
    echo  ""
    echo "  RUNNING  tRNA-Scan"
    #blastp -db /opt/PhageFinder/phage_finder_v2.5_4docker/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query "$1" -out ncbi.out -max_target_seqs 5 -num_threads 8 
    tRNAscan-SE --thread 1 -B -Q -o "${2}/tRNAscan.out" "$1"
    echo "  > tRNA-Scan DONE"
}

task_aragorn(){
    echo  ""
    echo "  RUNNING  Aragorn"
    /opt/Aragron/aragorn -m -o "${2}/tmRNA_aragorn.out" "$1";
    echo "  > Aragorn DONE"
}

task_phagefinder (){
    echo  ""
    echo "  RUNNING Phage_Finder_v2.5"
    /opt/PhageFinder/phage_finder_v2.5_4docker/bin/Phage_Finder_v2.5.pl -t "${1}/ncbi.out" -i "$2" -r "${1}/tRNAscan.out" -n "${1}/tmRNA_aragorn.out" -A "$3" -S
    echo "  > Phage_Finder_v2.5 DONE"
}

task_validationn (){
    if [ -s "${1}" ]; then
        echo "  VALIDATION: $1 exist and it not-empty"
    else
        echo "  VALIDATION: $1 doesn't exist or it is empty"
        exit 1
    fi
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

main() {
    echo  ""
    echo  "STARTING PHAGE FINDER SCRIPT"
    echo  ""
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

    # Validation 1
    [[ -z "${outputfile_location}" ]] && printf "    !!  Input Missing\n\n" >&2 && Help >&2 && exit 1

    # Echo inputs
    echo  "  INPUT faa_file: ${faa_file}"
    echo  "  INPUT fna_file: ${fna_file}"
    echo  "  INPUT phagefinder_info_file: ${phagefinder_info_file}"
    echo  "  INPUT outputfile_location: ${outputfile_location}"
    echo  ""

    # Validation 2
    echo  ""
    task_validationn "$faa_file" 
    task_validationn "$fna_file" 
    task_validationn "$phagefinder_info_file" 

    if [ ! -f "$outputfile_location" ]; then
        mkdir "$outputfile_location"
        echo "  DIR CREATED: $outputfile_location"
    fi
    

    # Run Tasks
    task_hmmsearch "${faa_file}" "${outputfile_location}"
    task_blast "${faa_file}" "${outputfile_location}"
    task_aragorn "${fna_file}" "${outputfile_location}"
    [[ "${transcan_file}" == "true" ]] && task_trnascan "${fna_file}" "${outputfile_location}"
    task_phagefinder "${outputfile_location}" "${phagefinder_info_file}" "${fna_file}"
    echo  ""
    echo  "PHAGE FINDER SCRIPT DONE"
    echo  ""
    
    return 0
}

main "${@:-}"
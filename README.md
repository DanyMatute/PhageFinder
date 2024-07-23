# PhageFinder

Phage_Finder was developed at The Institute for Genomic Research (TIGR) by Dr. Derrick E. Fouts as a heuristic computer program written in PERL to identify prophage regions within bacterial genomes (In press Nucleic Acids Research, 2006).

It uses tab-delimited results from NCBI BLASTALL or WU BLASTP 2.0 searches against a collection of bacteriophage protein sequences and results from HMMSEARCH analysis of 441 phage-specific HMMs to locate prophage regions. By using FASTA33, MUMMER or BLASTN, it can find potential attachment (att) sites of the phage region(s). Data from tRNAscan-SE and Aragorn are used to determine whether a tRNA or tmRNA served as the putative target for integration. Additionally, by looking for the presence or absence of specific proteins using specific HMM models, Phage_Finder can predict whether the region is most likely prophage and which type (Mu, P2, or retron R73), an integrated element, a plasmid, or a degenerate phage region.

The goal of this project is to provide an open-sourced, standardized and automated system to identify and classify prophages within prokaryotic genomes. It is hoped that this package will facilitate future studies on the biology and evolution of these prophages by providing a level of microbial genome annotation that was previously void.

Fouts DE. Phage_Finder: automated identification and classification of prophage regions in complete bacterial genome sequences. Nucleic Acids Res. 2006;34(20):5839-51. doi: 10.1093/nar/gkl732. Epub 2006 Oct 24. PMID: 17062630; PMCID: PMC1635311.



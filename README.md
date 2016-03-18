# check_miseq
Checks the composition of fastq files, typically from files generated by MiSeq runs

# checkMiSeq.py
- uses bowtie2 to map reads to a set of reference genomes
- outputs a csv formatted file specifiying the number of reads mapping to each reference and the number of unmapped reads for each sample
- input can be the path to a folder containing fastq files or a pathlist, which is a file containing a list of folders
- by default, the bowtie2 index of the reference genomes points to /home/veratai/RefSeqs/checkMiSeq_refseqs
- other defaults:
  minlen (minimum match length, CIGAR M) = 100, 
  minq (minimum mapping quality) = 0, 
  mins (minimum alignment score) = 0

- eg. command: 
python ~/scripts/checkMiSeq.py -path /home/veratai/test/ -output checkMiseq.txt -log checkMiseq.log -x /home/veratai/RefSeqs/checkMiSeq_refseqs

# checkMiSeq_refseqs_withoutHIVHCV.fasta
- fasta file of reference genomes
- from the initial version, the HIV and HCV genomes are removed - in order to integrate with the MiCall pipeline
- a bowtie2 index is required for the reference genomes to be used with checkMiSeq.py
- but it is not built, because these references will need to be combined with the mapping step in the MiCall pipeline
- the reference genomes (and their genbank accession) are:\n
  phiX174_sensulato_NC_001422, 
  EcoliK12_substrMG1655_NC_000913, 
  EcoliRR1_CP011113, 
  Pacnes_CP006032, 
  GBvirusC_NC_001710, 
  GBvirusB_NC_001655, 
  Human_pegivirus2_NC_027998, 
  HBV_ayw_NC_003977, 
  Moloney_murine_leukemia_virus_NC_001501, 
  Pdenitrificans_PD1222, 
  hg38 (human genome, version 38)

# checkMiSeq_allSamples.R
- R script to make barplots of the number of reads hitting the references from the csv file produced by checkMiSeq.py
- eg. command: 
Rscript ~/RScripts/checkMiSeq_allSamples.R checkMiseq.txt


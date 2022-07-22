# ZosteraSNPs

This is the SNP dataset for Zostera marina generated January 2020  

### map of sites and geographic distances
***R/SiteMaps.R***  
***data/NSDE41661/***  
***output/GeographicDistanceIndiv.csv***  
***output/IndivLatLon.csv***  
***output/SiteQuadratLatLon.csv***  
***output/SiteMaps_coastline.pdf***  
***output/SiteMaps.pdf***  

### map of sites and admix results
***R/SiteMaps_NGSadmix.R***  
***data/NSDE41661/***  
***output/SiteMaps_NGSadmix.pdf***  

### angsd-generated genotype likelihoods in beagle format  
***data/ind393_clean***  
***data/loci19433***  
***data/zos.393ind.HWE.99.gl.beagle.gz***  

### genotype calls
***R/SubsetGenotypeCalls.R***  
***data/zos.393ind.HWE.99.loci19433.vcf.gz***  
***data/genotypeCalls.393ind.8684loci.geno***This is the output dataset  
***data/ind393_clean***  
***data/loci8684***    

### Analysis of clones (as determined by ngsRelate) ###
***R/ngsRelate_ClonesOnMap.R***  
***R/Relatedness~Distance-ngsRelate.R***  
***data/ngsRelate/FILES***  The output from ngsRelate  
***data/ind245adults_noClones_clean***  The output of unique genets to use  
***ListOfClones.txt***  
***ngsRelate_ClonesOnMap=threshold=0.8.pdf***  
***output/Relatedness~Distance-ngsRelate-final.pdf***  

### Basic pop gen statistics ###    
***R/basicStats_adultsNoClones.R***  
***output/BasicStats_adultsNoClones.csv***  

### Create FASTA from calls and generate Maximum Likelihood tree ###  
***R/ConvertVCFToFASTA.R***  
***R/iQtreePlot.R***  
***data/fasta/XXX***  
***data/fasta/XXX***  
***data/fasta/XXX***  

### Analysis of Molecular Variance (AMOVA) ###
***R/amova_noClones_calls.R***  
***output/amova.method1-calls_noClones.csv***  
***output/amova.method1-calls_noClones-pairwisePhiSt.csv***  

### Admixture analysis - adults and seeds ###  
***R/ngsAdmix_dips.pretty_adults-NoClones.R***  
***R/ngsAdmix_dips.pretty_seeds.R***  
***data/ind245adults_noClones_clean***  
***data/ind99_seeds***  
***data/ngsadmix-noClones/*opt***  
***data/NGSadmix-seeds/*opt***  
***output/ngsAdmix_dips.pretty-adults-NoClones.pdf***  
***output/ngsAdmix_dips.pretty-seeds.pdf***  

# Commands from other programs #
### ngsRelate ###
***see https://github.com/ANGSD/NgsRelate***  
First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).***  
***angsd -b DD_43_bamfiles -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -rf zos.393ind.HWE.99.forANGSD.loci19433.txt -out DD***  
Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)  
***zcat DD.mafs.gz | cut -f5 |sed 1d > DD_freq***  
run NgsRelate The output should be a file called newres that contains the output for all pairs  
***./ngsRelate -g DD.glf.gz -n 43 -f DD_freq -O DD_newres***  

### bwa (Version: 0.7.17) and samtools (Version: 1.9) to align reads and create bamfiles ###
***### pull in lines from fastq list and perform bwa***  
***#! /bin/bash***  
***cat fastqfiles | while read LINE***  
***do***  
***bwa mem ../genome/zostera_genome_forBWA -t 20 ../fastq/$LINE.fastq | samtools sort - | samtools view -h - | samtools view -b > $LINE.bwa.bam***  
***samtools index $LINE.bwa.bam***  
***done***  

### bcftools for calling SNPs (Version: 1.9-259) ###
***bcftools mpileup -d 500 -Ou -a AD -f ../genome/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna -b ind393 |\***  
***bcftools call -mv -Ou -G Sample_group.txt |\***  
***bcftools view -m2 -M2 -q 0.01:minor -o zos.393ind.HWE.99.vcf***  
***bgzip zos.393ind.HWE.99.vcf***   
***bcftools index zos.393ind.HWE.99.vcf.gz***  
***bcftools query -f '%CHROM\t%POS\t%POS\n' zos.393ind.HWE.99.vcf.gz > zos.393ind.HWE.99.allloci.txt***  

### creating angsd-formatted beagle file (version: 0.931) ###
***angsd -GL 1 -out zos.393ind.HWE.99.geno -doGlf 2 -rf zos.393ind.HWE.99.allloci.forANGSD.txt -doMaf 2 -doMajorMinor 1 -SNP_pval 0.000001 -doGeno 5 -doPost 1 -bam ind393***  

### iQTree (version 1.6.11) ###
***Command: /usr/local/iqtree/bin/iqtree -s /tmp/IQTREE/_aGIwYTklJ/seqin.fst -st DNA -pre /tmp/IQTREE/_aGIwYTklJ/iqt -m GTR+F -bb 1000 -pers 0.5 -numstop 100 -nt 3***  

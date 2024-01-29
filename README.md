# ZosteraSNPs

This is the plant trait and SNP dataset for Zostera marina
Updated Dec 30 2023

### map of sites and geographic distances  
***R/SiteMaps.R***  
***data/NSDE41661/***  
***output/GeographicDistanceIndiv.csv***  
***output/IndivLatLon.csv***  
***output/SiteQuadratLatLon.csv***  
***output/SiteMaps_coastline.pdf***  
***output/SiteMaps.pdf***  

### Adult plant traits 
***R/DensityHeightAnalysis.R***  
***data/perm.quad.density.2019.csv***  
***data/perm.quad.veg.ht.2019.csv***  
***output/DensityHeight.png***  

### Genotype likelihoods and calls are in data/SNPs/###

### Analysis of clonality via ngsRelate
***data/ngsRelateFiles/***  The output from ngsRelate  
***R/ngsRelate_processFiles_rab.R***  
***output/histogram_ngsRelate_rab.pdf***  
***output/clones_r=70perc.txt***  
***output/clones_r=80perc.txt***  

### Analysis of clonality via dissimilarity  
***data/DissimilarityDist_clones_95perc.csv.gz***  
***R/dissimilarity.R***  
***output/histogram_dissimilarity.pdf***  

### finalize list of clones and make map ###  
***R/IdentifyClones_summary.R***  
***output/ListOfClones.txt***  
***output/Clone_analysis.csv***  
***data/214adults99seeds_noClones_clean***  The output of unique genets to use  
***R/ngsRelate_ClonesOnMap.R***  
***output/output/ClonesOnMap_Both.pdf***  

### Seed density and size   
***data/seedPhenotypes.csv***  
***R/seedPhenotypes.R***  
***output/seedPhenotypes.png***  

### map of sites and admix results
***R/SiteMaps_NGSadmix.R***  
***R/ngsAdmix_dips.pretty_adults-NoClones.R***  
***data/NSDE41661/***  
***data/ngsAdmixFiles/***  
***output/SiteMaps_NGSadmix.pdf***  
***output/ngsAdmix_dips.pretty-adults+seeds.pdf***  

### Analysis of Molecular Variance (AMOVA) and Pairwise PhiST
***R/amova_noClones_calls.R***  
***R/amova_seeds.R***  
***output/amova.method1-calls_noClones.csv***  
***output/amova.method1-calls_noClones-pairwisePhiSt.csv***  
***output/amova.method1-calls_seeds-pairwisePhiSt.csv***  
***output/amova.method1-seeds.csv***  

### Isolation by distance
***R/Fst~Distance.R***  
***output/amova.method1-calls_noClones-pairwisePhiSt.csv***  
***output/Fst~Distance.pdf***  

### Region-wide PCA ###  
***R/PCA_genotypeCalls_90perc.R***  
***output/PCA_genotypeCalls_90perc.pdf***  

### DAPC assignment - seeds & adults  
***R/SeedAssignment-DAPC.R***  
***output/SeedAssignment-DAPC.txt***  
***output/SeedAssignment-DAPC.pdf***  

### SPAGEDI correlogram of adults
***R/SpagediSGS_noClones_10categories.R***  
***data/data/allSite.Spagedi.out.txt***  
***output/SpagediSGS_noClones_10categories.pdf***  


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

# ZosteraSNPs

This is the SNP dataset for Zostera marina generated January 2020  

### angsd-generated genotype likelihoods  
***data/ind393_clean***  
***data/loci19433***  
***data/zos.393ind.HWE.99.gl.beagle.gz***  

### genotype calls pulled from vcf file
***R/SubsetGenotypeCalls.R***  
***data/zos.393ind.HWE.99.loci19433.vcf.gz***  
***data/genotypeCalls.393ind.8684loci.geno***  
***data/ind393_clean***  
***data/loci8684***  

### adult clone IDs ###
***data/ind245adults_noClones_clean***  
***data/CloneIDs.csv***   

### Heterozygosity ###    
***R/basicStats_adultsNoClones.R***  
***output/BasicStats_adultsNoClones.csv***  

### Create FASTA from calls and generate pretty ML tree ###  
***R/ConvertVCFToFASTA.R***  
***R/iQtreePlot.R***  
***data/fasta/393ind8683loc.fas***  
***data/fasta/245ind8683loc_noClones_foriQtree.fas***  
***data/fasta/99seeds_8683loc_foriQtree.fas***  



# Scripts for 
## Wang X, Pedersen CE, Athanasiadis G, Garcia-Erill G, Hanghøj K, Bertola LD, Rasmussen MS, Schubert M, Liu X, Li Z, Lin L. Persistent gene flow suggests an absence of reproductive isolation in an African antelope speciation model. bioRxiv. 2022:2022-12.

### 1.Mapping
A pipeline designed for mapping and post-mapping filtering of waterbuck.

### 2.Sample filtering 
A pipeline designed to remove problematic samples prior to downstream analyses, based on
`mapping statistics;`
-statistics from MultipleQC;
-exteremely high heterozygosity,
-error rates,
-and between-sample relatedness.

### 3.Sites filtering
A pipeline designed to avoid biases from low-quality mapping, based on 
detection of problematic regions of the reference genomes (mappability, repeats and sex-linked chromosomes or scaffolds);
sites that showed unusual depth; 
sites excess heterozygosity after mapping. 

### 4.Genotypelikehood calculation
Script to calculate genotype likelihood

### 5.Population structure and recently admixture
Script to infer PCA by PCAngsd;
Script to infer and plot admixture by NGSadmix;
Script to use evalAdmix;
Script to generate NJ-tree.

### 6.Ancient Admixture
Script to perform TREEMIX;
Script to calculate ABBABABA tests by ANGSD;

### 7.Divereisty and divergence
Script to infer Heterozygosity;
Script to infer Pairwise global Fst between each pair of populations;
Script to run EEMS;
Script to estimate divergence time and their 95%CI range between two-subspecies using Fastsimcoal27.

### 8.Genome-wide scans
Script to perform sliding window Fst;
Script to perform sliding window Dxy;
Script to perform sliding window TWISST;
Script to calculate sliding window neutralisty statistics;
Script to infer linkage disequilibrium by ngsLD;
Script to run EEMS;

### 9.Detection of gene flow within differentiation islands
Scripts for infer local population structure, calculate heterozygosity for each cluster and highlight wrong clustered samples in ABBABABA.

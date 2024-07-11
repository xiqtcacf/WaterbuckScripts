# Scripts for 
## Wang, X., Pedersen, C. E. T., Athanasiadis, G., Garcia-Erill, G., Hanghøj, K., Bertola, L. D., ... & Heller, R. (2022). Persistent gene flow suggests an absence of reproductive isolation in an African antelope speciation model. bioRxiv, 2022-12. https://www.biorxiv.org/content/10.1101/2022.12.08.519574v1.abstract
### 1.Mapping
A pipeline designed for mapping and post-mapping filtering of waterbuck.

### 2.Sample filtering 
A pipeline designed to remove problematic samples prior to downstream analyses, based on

mapping statistics;

statistics from MultipleQC;

exteremely high heterozygosity;

error rates;

relateness between pairwise samples: King calculated from global 2D-sfs;

relateness within each location by NgsRelate.


### 3.Sites filtering
A pipeline designed to avoid biases from low-quality mapping, based on 

3.1.detection of problematic regions of the reference genomes:

3.1.1.mappability;

3.1.2.repeats;

3.1.3.sex-linked chromosomes or scaffolds;

3.2.sites that showed unusual depth; 

3.3.sites excess heterozygosity after mapping. 


### 4.Genotypelikehood calculation
Script to calculate genotype likelihood.

### 5.Population structure and recently admixture
5.1.Script to infer Principal component analysis using PCAngsd;

5.2.Script to infer and plot admixture by NGSadmix;

5.3.Script to use evalAdmix;

5.4.Script to evaluate cases of recent admixture;

5.5.Script to generate NJ-tree.


### 6.Ancient Admixture
6.1.Script to perform TREEMIX;

6.2.Script to calculate ABBABABA tests by ANGSD;


### 7.Divereisty and divergence
7.1.Script to infer Heterozygosity;

7.2.Script to infer Pairwise global Fst between each pair of populations;

7.3.Script to run EEMS;

7.4.Script to estimate divergence time between two-subspecies using Fastsimcoal27;

7.5.Script to estimate divergence time's 95%CI range between two-subspecies using Fastsimcoal27;

7.6.All three demographic model design for estimating diveregence time by Fastsimcoal27;

7.7.Script to estimate divergence time using the Two-Two (TT) method.


### 8.Genome-wide scans
8.1.Script to perform sliding window Fst;

8.2.Script to perform sliding window Dxy;

8.3.Script to perform sliding window TWISST;

8.4.Script to calculate sliding window of nucleotide diversity and Tajima’s D;

8.5.Script to infer linkage disequilibrium by ngsLD;


### 9.Detection of gene flow within differentiation islands
Scripts for infer local population structure, calculate heterozygosity for each cluster and highlight wrong clustered samples in ABBABABA.

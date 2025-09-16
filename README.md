# Analysis of genetics underlying the hybrid effect variation in a Eurasian pig cross
Analysis carried out by Hai-Bing Xie (xiehb@mail.kiz.ac.cn), Zi-Qin Huang, Li-Gang Wang, and Long-Chao Zhang in consultation with Ya-Ping Zhang, Chung-I Wu and Li-Xian Wang.

# Overview

This provides the data and scripts to explore the genetics underlying hybrid effect variation in a F2 population of the pig (*Sus scrofa*). The F2 was developed from a cross between the European Large White (LW) boars and East Asian Min (MIN) sows.

## 1. The phenotypic data for the LW-MIN F2 population
All the F2 individuals were raised without any directional selection imposed by humans, and the phenotypes of the F2 are expected to vary with the segregation of the LW and MIN alleles in the F2 population. A total of 135 traits were collected in the F2 population, and most of them were collected at 240 days after birth. The raw phenotypic data was supplied in CSV format with a filename of “f2_trait_name_trait_value.csv”. The file contains five columns, indicating the sample id ("f2"), family id ("family"), an assigned number of a trait ("trait_id"), the abbreviation of the trait name ("trait_name"), and the phenotypic value of the trait ("trait_value"). The male sex is given with odd numbers of sample id and female sex with even numbers. The following shows an example of the intramuscular fat content (IMF) trait for a subset of the F2 samples.

```
f2	family	trait_id	trait_name	trait_value
1008001	10080	67	IMF	1.53
1008905	10089	67	IMF	1.87
1008902	10089	67	IMF	1.13
1008003	10080	67	IMF	3.13
932304	9323	67	IMF	1.4
932302	9323	67	IMF	2.47
932403	9324	67	IMF	2.81
932412	9324	67	IMF	2.9
931809	9318	67	IMF	2.16
932307	9323	67	IMF	6.31
```

## 2. Genotyping and phasing
For all the F0, F1, and F2, the genomes were determined using the Illumina SNP60 BeadChip. The genomic data is available in the file "sscrofa11.data.zip". The coordinates were shown in the genome build of Sscrofa 11.1.

We performed the haplotype inference for the whole LW-MIN family using the  shapeit software.

```
for i in `seq 1 18`
do
	shapeit -T 10 -P chr$i.ped chr$i.map --duohmm -W 15 --output-max chr$i.phased.duohmm --output-graph chr$i.phased.duohmm.graph --force
done
```
* [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) - Version : v2.r904


## 3. Transmission analysis and recombination breakpoint identification
We prepared a R script ("1.determination of recombination breakpoints and allelic transmission.R") to determine the genomic transmission from F1 to F2 and the inheritance of LW or MIN alleles in the paternal and maternal genomes of each F2. 

We applied several steps to determine the genomic transmission in the LW-MIN family. First, the phased data produced by the shapeit software were loaded into a dataframe in  the R. To reduce the computation complexity, the haplotype data of four grandparents and two parents for each F2 were extracted for analysis. Second, a combination of R functions *GetScore*, *StepHap*, and *correcting_haplotype* were developed to determine the genomic transmission across generations, especially from F1 to F2. In the LW-MIN family,  the F1 is always heterozygous (paternal/maternal as LW/MIN). In F1 males, the paternal chromosome (LW) was assigned with a number of *5*, and maternal chromosome (MIN) with a number of *6*. The numbers do not indiciate the chromosomes 5 and 6, but show the allele origin in the F0-F1-F2 transmission for each autosome. In F1 females, the paternal and maternal chromosomes were assigned to *11* and *12*. For the F2 offspring, the paternal and maternal chromosomes  were assigned *3* and *4*, respectively. Therefore, the *3* was identical to *5* or *6* if no recombination events was involved in the paternal genome, or a mosaic chromosome with both *5* and *6*. A similar case is for *4* represented by the *11* and *12*.  The following figure shows the naming convention.

![image](https://github.com/xiehb-evolution/hybrid-effects/blob/main/tmp/inheritance2.jpg)


The  *5*, *6*, *11*, and *12* were designed to indicate the inheritance of LW (odd number) or MIN (even number) alleles in the paternal and maternal genomes of the F2. The recombination events were determined by finding a boundary haborbing allele changes, for example from *5* to *6* on the paternal chromosome with an assigned number *3*. 

For details, please see the implementation of the three R functions.

The genomic transmission for all SNPs were generated in the f2.inheritance.txt file with the following format:
```
f2	chromosome	position	paternal_allele	maternal_allele
1007207	2	137601377	5	11
1007207	2	137609824	5	11
1007207	2	137650460	5	11
1007207	2	137892857	5	11
1007207	2	137978735	5	12
1007207	2	138011068	5	12
1007207	2	138014980	5	12
1007207	2	138048646	5	12
1007207	2	138066066	5	12
1007207	2	138137424	5	12
```
There is a recombination event detected in the maternal chromosome 2 of the F2 (id 1007207) between two SNPs (coordinates:137892857 and 137978735). The fifth and sixth columns indicates the F0 alleles in the F2 genomes as described above (*5*: LW allele, *11*: LW allele, and *12*: MIN allele). 

The following figure shows the inheritance of alleles on chromosome 5 in a F2 individual (930806). The red and blue colors indicate the LW and MIN alleles, respectively. The recombination breakpoints are indicated by the boundary SNPs at red-blue color shift  in the paternal and maternal genomes.

![image](https://github.com/xiehb-evolution/hybrid-effects/blob/main/tmp/shapeit.jpg)

To reduce the complexity in the following analysis on hybrid effects, the f2.inheritance.txt was saved to a MySQL server, all the subsequent analysis was based on the saved tables.

## 4. Hybrid effect analysis

### 4.1 Data processing and window-based analysis

The inheritance data from Section 3 was processed using the R script "2. F2 inheritance fragment processing and trait association preparation.R" to convert SNP-level inheritance patterns into 100-kb genomic windows. For each window, F2 individuals were classified into four genotype categories:

- **MIN/MIN**: Both paternal and maternal alleles from MIN population
- **LW/LW**: Both paternal and maternal alleles from LW population
- **MIN/LW**: Paternal MIN allele, maternal LW allele
- **LW/MIN**: Paternal LW allele, maternal MIN allele

The analysis maintains distinction between reciprocal heterozygotes to account for parent-of-origin effects.

### 4.2 Definition of hybrid effects

The hybrid effect was defined as the difference in phenotypic means of the homozygous and heterozygous genotypes in the LW-MIN F2 population. Within each 100-kb window, phenotypic means were calculated for each genotype on the 135 traits. The analysis was conducted separately for males and females. Four homozygote-heterozygote comparisons were performed:

1. MIN/MIN vs MIN/LW
2. MIN/MIN vs LW/MIN
3. LW/LW vs MIN/LW
4. LW/LW vs LW/MIN

The hybrid effect size was calculated as below:

```
Standardized hybrid effect size = abs(Mean-of-heterozygote - Mean-of-homozygote) / (Mean of a sex)
```

This convertion enables direct comparison across all 135 traits with different units and scales.

### 4.3 Exponential distribution of hybrid effect size and λ parameter estimation

The distribution of standardized hybrid effect sizes was fitted to exponential distributions using the `fitdistrplus` R package. 

To accurately investigate the hybrid effect in sexes, the λ parameter was estimated separately for male and female individuals. The hybrid effect size distributions are presented in the figure below, with figure A indicating the distribution in males and figure B indicating the distribution in females.

![image](https://github.com/xiehb-evolution/hybrid-effects/blob/main/tmp/figureS1_00.png)

### 4.4 Interplay of genetic differentiation and recombination in hybrid effect transition
The hybrid effect variation was explored across 100-kb autosomal windows under different levels of genetic differentiation (Fst). The genomic windows were classified into 20 bins according to their Fst values calculated between the LW and MIN founder populations using the whole genome resequencing data. The λ for the exponential distribution of hybrid effects and the mean recombination rate were calculated for each bin of the autosomal windows. The genomic transition of hybrid effects were examined according to the variation of the λ parameters across bins. The following figure shows a hump-shaped distribution of the λ and a positive correlation between recombination rates and λ. The R script is provided in the file "figure1.R".

![image](https://github.com/xiehb-evolution/hybrid-effects/blob/main/tmp/lambda.png)

According to Ronald Fisher's geometric model (FGM, please see **THE GENETICAL THEORY OF NATURAL SELECTION**, 1930), mutations with large effect sizes tend to be deleterious. In the FGM framework, we explored the hybrid effect by considering the homozygotes as wild types and the heterozygotes as mutant types. The outcome of hybridization is expected to be affected by the hybrid effect size defined in section **4.2**. The expected effect size for each of the 20 bins is the reciprocal value of the λ parameter (or 1/λ), and therefore larger λ values indicate higher probabilities of being beneficial hybrid effect in the FGM. 

Here, we chose an empricial threshold at λ = 33.5 (close to the λ mean for the 20 bins) to define three distinct types of hybrid effect:
- **Hybrid vigor**: λ > 33.5 (small effects in regions showing mild-to-moderate levels of differentiation and the highest recombiantion rates)
- **Inbreeding depression**: λ < 33.5 and Fst < 0.095 (large effects in weakly differentiated regions)
- **Hybrid depression**: λ < 33.5 and Fst > 0.285 (large effects in highly differentiated regions)


### 4.5 Sex difference in hybrid effect variation

The λ values were calculated independently for F2 males and females to reveal the sex difference in hybrid effect variation. The constitutively lower λ values across different bins pinpoint an evolutionary inferior position of females in hybridization. The following figure shows the sex difference in the distribution of λ parameters. The R script is provided in the file "figure2.R".

![image](https://github.com/xiehb-evolution/hybrid-effects/blob/main/tmp/figure2_A.jpg)


### 4.6 Database implementation and quality control

All data was integrated into MySQL databases using the following workflow:

1. **Data integration**: The "f2inheritance.txt" file from Section 3 was imported into MySQL tables
2. **Fragment processing**: SNP-level data was aggregated into 100-kb windows
3. **Trait association**: Phenotypic data was linked to genotype classifications
4. **Statistical computation**: Window-based means and effect sizes calculated for each trait-sex combination

Quality control filters excluded windows with insufficient sample sizes or excessive missing data. The database schema was optimized for efficient querying of the large-scale genomic and phenotypic datasets.

The key database tables created during the analysis include:

**f2inheritance**: SNP-level allelic inheritance patterns

```
f2	chr	pos	f1father	f1mother
1007207	2	137601377	5	11
1007207	2	137609824	5	11
1007207	2	137650460	5	11
```

**f2fragmentinheritance_window100k**: 100-kb window genotype assignments

```
f2	chr	window	inheritance	origin
1007207	1	0	0	P
1007207	1	0	1	M
1007207	1	1	1	P
```

**window100k_single_site_trait_stat_mutant_deviation_from_mean**: Hybrid effect statistics

```
trait_id	chr	window	paternalinheritance	maternalinheritance	malecount	femalecount	maletrait	femaletrait	maledev	femaledev
1	1	0	0	0	15	18	2.34	2.56	0.12	0.08
1	1	0	0	1	12	14	2.45	2.61	0.18	0.14
1	1	0	1	0	13	16	2.29	2.48	0.09	0.05
```

**renew_complete_data_with_overall**: Integrated analysis table with FST and recombination data

```
chr	window	lambda	WEIGHTED_FST	snps	recombination_rate	Sex	tag
1	0	35.2	0.124	45	0.0034	Overall	hybrid vigor
1	1	28.7	0.089	52	0.0041	Overall	inbreeding depression
1	2	31.8	0.267	38	0.0028	Overall	hybrid depression
```

The computational pipeline enables genome-wide analysis of hybrid effects while maintaining efficiency through the sliding window approach. All R scripts for data processing, statistical analysis, and visualization are provided in this repository.


## 5. Female heterozygote deficiency
The whole genome resequencing data was deposited in the Genome Sequence Archive (GSA; http://gsa.big.ac.cn) under accession number CRA002451. The lists for F2 males and F2 females are provided in the text files "F2male.txt" and "F2female.txt".

To reduce the detection of heterozygotes caused by unreliable/multiple mapping of paralogous reads in the autosomes, the BAM files (mapping data of the genomic reads) were filtered using a C++ program ("filter_bam.cpp" under the folder "Cpp") to keep uniquely-mapped paired-reads with high mapping qualities. The compilation of the executable binary program requires the HTSLib. The usage of filter_bam is supplied with an example of shell script "process_bam.sh".
* [HTSlib](https://github.com/samtools/htslib/releases/) - Version : v1.20

The comparison between male and female heterozygote frequencies is used to reveal the selection in each sex. The analysis is based on the whole genome resequencing data. The SNPs were grouped by the minior allele frequency (MAF) with a bin size of 0.05. For each group, the number of heterozogytes and homozygotes were compared between the sexes. 

The analysis was conducted using a C++ program (xie_unphased_vcf_for_heterozygote_stat.cpp under folder "Cpp"):
```
vcftools --gzvcf F2.biallelic.chr.vcf.gz --min-alleles 2 --max-alleles 2 --exclude-positions excluded.snps.list --recode --stdout | xie_unphased_vcf_for_heterozygote_stat 100000 4 15 F2male.txt F2female.txt > F2.4to15X.100k.exclude.snps.allelefreq.stat.out
```
Here, the 100000 is the size of sliding windows, 4 and 15 are the thresholds of minimum and maximum sequencing depths for SNP sites, and the F2male.txt and F2female.txt are files providing the sample lists of F2 males and F2 females. The "F2.biallelic.chr.vcf.gz" is the VCF file for the LW-MIN family. The "excluded.snps.list" is a list of SNPs that are excluded in the analysis due to the tendency of mapping errors from paralogous genomic sequences.

The results is provided in the "F2.4to15X.100k.exclude.snps.allelefreq.stat.out" file. The format is given as below:
```
Chr     Position        Pop1size        Pop2size        SNPs    HomoSites1      HetSites1       HomoSites2      HetSites2       HetRatio1       HetRatio2
```
The Chr and Position indicate a 100-kb sliding window with a given lower boundary (0-based), the Pop1size (from F1male.txt) and Pop2size (from F2female.txt) are the sizes of samples with genomic data covering this window. SNPs is the total number of SNPs in this window with 0 < MAF < 0.05. HomoSites1 and HetSites1 are the total number of homozygotes and heterozygotes counted in the samples in F2male.txt, and HomoSites1 and HetSites1 provide information in samples from F2female.txt. HetRatio1 and HetRatio2 are the ratio of heterozygotes in the F2 males and females in a window. The seven columns (from SNPs to HetRatio2) replicate 10 times to indicate the results on SNPs with different MAFs (a step size of 0.05).







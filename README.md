#Analysis of the genetics underlying the hybrid effect variation in a Eurasian pig cross
Analysis carried out by Hai-Bing Xie (xiehb@mail.kiz.ac.cn) in consultation with Ya-Ping Zhang, Chung-I Wu and Li-Xian Wang.

# Overview

This provides the data and scripts to explore the genetics underlying hybrid effect variation in a F2 population of the pig (*Sus scrofa*). The F2 was developed from a cross between the European Large White (LW) boars and East Asian Min (MIN) sows.

## 1. The phenotypic data for the LW-MIN F2 population
All the F2 individuals were raised without any directional selection imposed by humans, and the phenotypes of the F2 are expected to vary with the segregation of the LW and MIN alleles in the F2 population. A total of 135 traits were collected in the F2 population, and most of them were collected at 240 days after birth. The raw phenotypic data was supplied in CSV format with a filename of “f2_trait_name_trait_value.csv”. The file contains five columns, indicating the sample id ("f2"), family id ("family"), an assigned number of a trait ("trait_id"), the abbreviation of the traits name ("trait_name"), and the value of the trait ("trait_value"). The male sex is given with odd numbers of sample id and female sex with even numbers. The following show an example for the intramuscular fat content (IMF) trait for a subset of the F2 samples.

```
f2	family	trait_id	trait_name	trait_value
1008004	10080	67	IMF	1.07
1008908	10089	67	IMF	2.48
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
932306	9323	67	IMF	3.4
932408	9324	67	IMF	2.28
```

## 2. Genotyping and phasing
For all the F0, F1, and F2, the genomes were determined using the Illumina SNP60 BeadChip. The genomic data is available in the file "sscrofa11.data.zip".

We performed the haplotype inference for the whole LW-MIN family using the  shapeit software.

```
for i in `seq 1 18`
do
	shapeit -T 10 -P chr$i.ped chr$i.map --duohmm -W 15 --output-max chr$i.phased.duohmm --output-graph chr$i.phased.duohmm.graph --force
done
```

* [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) - Version : v2.r904


## 3. Transmission analysis and recombination breakpoint identification
We prepared a R script (with a filename of "1.determination of recombination breakpoints and allelic transmission.R") to determine the genomic transmission from F1 to F2 and the inheritance of LW or MIN alleles in the paternal and maternal genomes of each F2. The analysis was performed on 100-kb autosomal windows.

We applied several steps to determine the genomic transmission in the LW-MIN family. First, the phased data produced by the shapeit software were loaded into a dataframe in  the R. To reduce the computation complexity, the haplotype data of four grandparents and two parents for each F2 were extracted for analysis. Second, a combination of R functions *GetScore*, *StepHap*, and *correcting_haplotype* were prepared to determine the genomic transmission across generations, especially from F1 to F2. In the LW-MIN family,  the F1 is always heterozygous (paternal/maternal as LW/MIN). In F1 males, the paternal chromosome (LW) was assigned with a number of *5*, and maternal chromosome (MIN) in F1 with a number of *6*. The numbers do not indiciate the chromosomes 5 and 6, but show the allele origin in the F0-F1-F2 transmission. In F1 females, the paternal and maternal chromosomes were assigned to *11* and *12*. For the F2 offspring, the paternal and maternal chromosomes  were assigned *3* and *4*, respectively. Therefore, the *3* was identical to *5* or *6* if no recombination events was involved in the paternal genome, or a mosaic chromosome with both *5* and *6*. A similar case is for *4* represented by the *11* and *12*.

The  *5*, *6*, *11*, and *12* were designed to indicate the inheritance of LW (odd number) or MIN (even number) alleles in the paternal and maternal genomes of the F2. The recombination events were determined by finding a boundary haborbing allele changes, for example from *5* to *6* on the paternal chromosome (*3*). 

For details, please see implementation of the three R functions.

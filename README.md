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

## 2. Phasing
For all the F0, F1, and F2, the genomes were determined using the Illumina SNP60 BeadChip. The genomic data is available in the file "sscrofa11.data.zip".

We performed the haplotype inference for the whole LW-MIN family using the  shapeit software.

```
for i in `seq 1 18`
do
	shapeit -T 10 -P chr$i.ped chr$i.map --duohmm -W 15 --output-max chr$i.phased.duohmm --output-graph chr$i.phased.duohmm.graph --force
done
```

* [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) - Version : v2.r904 



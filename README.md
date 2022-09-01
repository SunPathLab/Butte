# Butte ![Butte][badge_butte]


```
                                           ___________
                                           |    7    | __
                                           |    6    | ||
            ________                       |    5    | ||
            |      |                       |    4    | ||
           |       |                       |CN  3    ~~ |
          /         \                     /     2        \ 
_________/           \___________________/      1         \__________
```


BoUnds of Time Till Expansion, a computational framework to investigate the SCNA arrival time in the somatic evolution toward the most recent common ancestor of tumor sample(s).


## Installation

In R environment, make sure the following dependencies are installed:
```sh
lpSolve
GenomicRanges
IRanges
S4Vectors
dplyr
```

then
```sh
library(devtools)
install_github("SunPathLab/Butte")
```

## Getting Started
One can analyze the timing by running the function `Butte()` . The following is an example usage of the function Butte:

First, load the example input data provided in /demo.
```sh
library("Butte")
load("exampleData.rda") #This example data is provided in the directory "demo"
```
exampleData.rda have two dataset: exampleData_CN41 and exampleData_CN62, which are
data of two different segments with copy number states 4:1 and 6:2 respectively. In this section, we will use exampleData_CN62.

Specify following values:
```sh
point_mutation_read_counts = exampleData_CN62$x #vector
read_depths = exampleData_CN62$m #vector
sample_purity = exampleData_CN62$purity
history_matrices = cnmutHistory(6,2) #Total copy: Minor copy = 6:2
timing_type = "butte" #since the copy number state is 6:2 which has multiple possible histories. For identifiable histories, the type of timing estimation should be "identifiable".
```
Now run the function Butte:
```sh
Result = Butte(x=point_mutation_read_counts , m=read_depths , history=history_matrices , nt=6, nb=2, qmethod="fullMLE", type = timing_type , bootstrapCI="bootstrap", purity=sample_purity , B=100)
```
One can get the timing estimation from `Result$pi`. Since in 6:2 is an unidentifiable CN state, an upper and lower bound of timing are hence provided instead. `Result$pi` provides the confidence interval of the estimated bounds.


## Batch Analysis
We also provides a batch version of SCNA timing analysis, `scnaTiming()`.

An example usage is shown below:
```sh
library(Butte)
#run scna Timing
butte.res = scnaTiming(scnaFile=scnaFile, ssnvFile=ssnvFile, sn=sampleName, outname=sampleName,
                                                   public=TRUE, skipchunk = 100, B=500, pubOrSub="pubOrSub")

#save results
saveRDS(butte.res, file=paste0(outdir, "/butte.res.", sampleName, ".rds"))
write.table(butte.res[[2]], file=paste0(outdir, "/butte.tab.", sampleName, ".tsv"), sep="\t", quote=F, row.names=F)
```

`scnaFile` and `ssnvFile` are the two input files, respectively. See example files in the folder demo.

`scnaFile` is a tab delimited file of SCNA (somatic copy number alteration) segments as rows. The columns are defined as follows:
```sh
chrom: chromosome id
loc.start: segment start (coordinate)
loc.end: segment end
num.mark: number of SNPs used for SCNA calling for the segment
seg.mean: the raw averaged read ratio between tumor and normal
copynumber: total copy number
minor_cn: copy number of the minor allele 
major_cn: copy number of the major allele
allelicratio: read count of the reference allele divided by total depth
LOHcall: interpretable TITAN state (see R package TitanCNA for details)
cellularprevalence: proportion of tumor cells having this SCNA
ploidy: average tumor ploidy (haploid coverage factor)
normalproportion: the fraction of normal contamination
logcopynumberratio: log2 transformed copy number ratio between tumor and normal
```

`ssnvFile` is a tab delimited file with SSNVs (somatic single nucleotide variants) as rows. The columns are defined as follows:
```sh
Assuming XXX is the name of a sample:
chr: chromosome id
pos: the coordinate of the SSNV
XXXmafc: allele frequency (alternative read count/total depth)  
XXXaltc: alternative read count  
XXXrefc: reference read count  
XXXpu: the tumor purity of the sample  
XXXccf: cancer cell fraction of the SSNV  
XXXccfSD: standard deviation of ccf  
pubOrSub: if the mutation is public or private (subclonal).
```


## Reference

Zicheng Wang, Yunong Xia, Lauren Mills, Athanasios N. Nikolakopoulos, Nicole Maeser, Jason M. Sheltzer and Ruping Sun (2022)
**Evolving copy number gains promote tumor expansion and bolster mutational diversification**
[bioRxiv, 2022.06.14.495959, Jun. 2022](https://doi.org/10.1101/2022.06.14.495959)


## Contact
Ruping Sun: ruping@umn.edu
Yunong Xia: xia00045@umn.edu
Zicheng Wang: wang2569@umn.edu

[badge_butte]:      assets/badges/badge_butte.svg

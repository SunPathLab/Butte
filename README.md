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

## Intro

Butte is a computational framework to calculate the arrival time and initiation time of a clonal somatic copy number alteration (SCNA) observed in patient tumor sample(s). Currently the method focuses on copy number gains. A genomic region at an observed clonal SCNA state has evolved during the timeline from the germline to the founder cell of the clonal expansion. The timeline can be divided into three fractions in term of the copy number evolution. The first time fraction (T0) is the SCNA initiation time when the first gain occurs. The third fraction (TK) is the arrival time which measure the delay from the last gain to the start of population expansion. Butte can estimate the initation time and arrival time of an SCNA (with the total copy number as high as 7) given the read counts of somatic single nucleotide variants (SSNVs) occurred within corresponding genomic region and tumor purity. To do so, Butte adopts EM algorithm to find the allele state distribution of SSNVs. Then it either directly solve the time fraction (for SCNAs with identifiable history matrices), or adopts linear programming to calculate the upper bounds of these time fractions if multiple history matrices can exist for the SCNA or the underlying linear system is underdetermined.

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
One can analyze the timing by running the function `Butte` . The following is an example usage of the function Butte:

First, load the example input data provided in /demo.
```sh
library("Butte")
load("exampleData.rda") #This example data is provided in the directory "demo"
```
exampleData.rda has two list objects: `exampleData_CN41` and `exampleData_CN62`, which contains the input data for two SCNA segments with copy number states at `4:1` and `6:2`, respectively. For example,

Specify following values:
```ruby
totalCN = 6
minorCN = 2
example = exampleData_CN62
mutant_read_depth = example$x
total_depth = example$m
tumor_purity = example$purity
type = "butte"    #type should be set to "identifiable" for the following 

Result = Butte(x=mutant_read_depth, m=total_depth, history=history_matrices, nt=totalCN, nb=minorCN, qmethod="fullMLE", bootstrapCI="bootstrap", purity=tumor_purity, B=100)
```

SCNA `6:2` has multiple possible histories, `Butte` calculates the upper bounds of estimated timing. For `4:1` (in `exampleData_CN41`), the timing will be solved exactly instead. `Result$pi` outputs the upper bounds of the arrival time (pK) and initiation time (p0). `Result$piCI` provides the confidence interval of the upper bounds.


## Batch Analysis
We also provides a batch version of SCNA timing analysis, `scnaTiming()`.

An example usage is shown below:
```ruby
library(Butte)
#run scnaTiming

butte.res = scnaTiming(scnaFile=scnaFile, ssnvFile=ssnvFile, sn=sampleName, outname=sampleName,
                                                   public=TRUE, skipchunk = 100, B=100, pubOrSub="pubOrSub")

#save results
saveRDS(butte.res, file=paste0(outdir, "/butte.res.", sampleName, ".rds"))
write.table(butte.res[[2]], file=paste0(outdir, "/butte.tab.", sampleName, ".tsv"), sep="\t", quote=F, row.names=F)
```

`scnaFile` and `ssnvFile` are the two input files, respectively. See example files in the folder demo.

`scnaFile` is a tab delimited file of SCNA (somatic copy number alteration) segments as rows. The columns are defined as follows:
```
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
```
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

see [the pdf manual of functions in Butte](https://github.com/SunPathLab/Butte/blob/main/man/Butte_0.0.0.9000.pdf) 


## Reference

Zicheng Wang, Yunong Xia, Lauren Mills, Athanasios N. Nikolakopoulos, Nicole Maeser, Jason M. Sheltzer and Ruping Sun (2022)
**Evolving copy number gains promote tumor expansion and bolster mutational diversification**
[bioRxiv, 2022.06.14.495959, Jun. 2022](https://doi.org/10.1101/2022.06.14.495959)


## Contact
Ruping Sun: ruping@umn.edu
Yunong Xia: xia00045@umn.edu
Zicheng Wang: wang2569@umn.edu

[badge_butte]:      assets/badges/badge_butte.svg

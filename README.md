

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
Given a genome segment, one can analyze the timing by running the function `Butte()` on tumor data. The following is an example usage of function Butte:

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
`scnaFile` and `ssnvFile` are two input files of SCNA and SSNV. Example files are provided. For users who want to try the batch functionality on their own data, please convert data to the format of the example data. 

`scnaFile` is a table with segments as rows and other information as columns. The columns are defined as follows:
```sh
chrom: chromosome of the segment
loc.start: starting point of the segment according to hg38
loc.end: end point of the segment according to hg38
num.mark: number of somatic point mutations on the segment
seg.mean:
copynumber: total copy number
minor_cn: copy number of minor allele 
major_cn: copy number of major allele
allelicratio:
LOHcall:
cellularprevalence: proportion of cells out of the entire sample that contains the event
ploidy: 
normalproportion: 
logcopynumberratio: 
```
`ssnvFile` is a table with point mutations as rows and other information as columns. The columns are defined as follows:
```sh
Assuming XXX is the name of a sample:  
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
[Paperlink](link to be added)


## Contact
Ruping Sun: ruping@umn.edu
Yunong Xia: xia00045@umn.edu
Zicheng Wang: wang2569@umn.edu

[badge_butte]:      assets/badges/badge_butte.svg

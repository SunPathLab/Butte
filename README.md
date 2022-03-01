# Butte

```
                                         ___________
                                         |    7    |
                                         |    6    |
                                         |    5    |
                                         |    4    |
                                         |CN  3    |
                                        /     2     \ 
_______________________________________/      1      \__________
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

Example:
```sh
library(Butte)

#run scna Timing
butte.res = scnaTiming(scnaFile=scnaFile, ssnvFile=ssnvFile, sn=sampleName, outname=sampleName,
                                                   public=TRUE, skipchunk = 100, B=500, pubOrSub="pubOrSub")

#save results
saveRDS(butte.res, file=paste0(outdir, "/butte.res.", sampleName, ".rds"))
write.table(butte.res[[2]], file=paste0(outdir, "/butte.tab.", sampleName, ".tsv"), sep="\t", quote=F, row.names=F)
```

Format of input files:
```sh
```

## Reference

Zicheng Wang, Yunong Xia, Lauren Mills, Athanasios N. Nikolakopoulos, Nicole Maeser, Jason M. Sheltzer and Ruping Sun (2022)
*Evolving copy number gains promote tumor expansion and bolster mutational diversification*
[Paperlink](link to be added)



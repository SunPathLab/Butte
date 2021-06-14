#' estimate Timing for each SCNA segment 
#' 
#'
#' @param scnaFile the SCNA segmentation file
#' @param ssnvFile the SSNV file
#' @param sn sample name
#' @param outname output timing file name (including path)
#' @param public somatic timeline is focused on public mutations (across multi-samples) or not (otherwise, clonal variants in the specific sample)
#' @param pubOrSub the colname for column indicating if the mutation is public or not
#' @param mmut minimum number of mutations for running timing analysis
#' @param qmethod the method for estimating q (probabilities of a randomly acquired mutation having allele state of aj/Nt)
#' @param skipchunk segments with number of data points (probes) no more than this number will be skipped 
#' @return list: timing result; timing table (for visualization) and merged CNA data frame
#' @export
scnaTiming <- function(scnaFile, ssnvFile, sn, outname, public=FALSE, pubOrSub="pubOrSub",
                      skipchunk = 19, mmut=10, qmethod="fullMLE") {

    timingInput = cnmutData(scnaFile=scnaFile, ssnvFile=ssnvFile, skipchunk=skipchunk)
    cnvA2   = timingInput[[1]]
    sampAB  = timingInput[[2]]
    cnvHits = timingInput[[3]]
    snvHits = timingInput[[4]]
    
    #produce results
    result = list()
    resultTable = vector()
    li = 1
    for (i in 1:dim(cnvA2)[1]) {
        cnchrom = as.character(cnvA2$chrom[i])
        cnstart = cnvA2$loc.start[i]
        cnend = cnvA2$loc.end[i]
        cnLOHcall = as.character(cnvA2$LOHcall[i])
        cnminor = cnvA2$minor_cn[i]
        cnmajor = cnvA2$major_cn[i]
        cntotal = cnvA2$copynumber[i]
        message(paste(cnchrom, cnstart, cnend, cntotal, cnminor, sep="\t"))

        #skip SCNAs with Nt >= 8
        if (cntotal >= 8 | (cntotal == 2 & cnminor == 1)) {
            message(paste("skip cnv", i, "Nt >= 8 or normal diploid", dim(chopped)[1], sep=" "))
            next
        }
        
        # chop data to study the relevant ssnvs inside the SCNA region
        chopped = sampAB[snvHits[which(cnvHits == i)],]
        if (public == TRUE) {    #focus on the truncal ones to multi samples
            chopped = chopped[which(chopped[,pubOrSub] == "public"),]                                    
        } else {                 #focus on the clonal ssnvs in the corresponding sample
            chopped = chopped[which(chopped[,paste(sn,"ccf",sep="")] +
                                    2.58*chopped[,paste(sn,"ccfSD",sep="")] >= 1),]
        }

        #skip SCNAs containing less than mmut mutations
        if (dim(chopped)[1] < mmut) {
            message(paste("skip cnv", i, "< min # of mutations", dim(chopped)[1], sep=" "))
            next
        }
        
        purity = chopped[,paste(sn,"pu",sep="")][1]
        
        onlyMuts = data.frame(chromosome=chopped$chr, position=chopped$pos,
                              refbase=chopped$ref, mutbase=chopped$alt, rsID=chopped$id,
                              t_ref_count=chopped[,paste(sn,"refc",sep="")],
                              t_alt_count=chopped[,paste(sn,"altc",sep="")],
                              allelefreq=chopped[,paste(sn,"mafc",sep="")],
                              t_depth=chopped[,paste(sn,"refc",sep="")] +
                                  chopped[,paste(sn,"altc",sep="")])
                              #t_depth=chopped[,paste(sn,"d",sep="")])
        
        # SCNA has two groups
        # group1: one group with an unique and identifiable history matrix A
        hmatrix = cnmutHistory(nt = cntotal, nb = cnminor)
        if (length(hmatrix) == 1 & dim(hmatrix[[1]])[1] == dim(hmatrix[[1]])[2])  #2:0 3:1 3:0 4:1
            cntype = "identifiable"
        
        # group2: the other does not have identifiable history, or there are mulitple A
        else
            cntype = "butte"
        message(cntype)

        # call Butte for estimating timing
        x <- Butte(x=onlyMuts$t_alt_count, m=onlyMuts$t_depth, history=hmatrix, qmethod=qmethod,
                   nt=cntotal, nb=cnminor, minMutations=mmut, type=cntype,
                   purity=purity, bootstrapCI="bootstrap", B=100)
        
        x = c(x, cnid=i, cnchrom=cnchrom, cnstart=cnstart, cnend=cnend,
              cnmajor=cnmajor, cnminor=cnminor, cnLOHcall=cnLOHcall)

        result[[li]] = x
        names(result)[[li]] = i

        K = length(x$pi)
        currentline = c(cnchrom, cnstart, cnend, cnmajor, cnminor, x$summaryTable[1], cntype,
                        as.numeric(x$pi[1]), x$piCI[1,1], x$piCI[1,2],
                        as.numeric(x$pi[K]), x$piCI[K,1], x$piCI[K,2])
        names(currentline) = c("chrom", "loc.start", "loc.end", "major_cn", "minor_cn", "nmut", "type",
                               "p0","p0l","p0h","pK","pKl","pKh")
        
        message(paste(currentline, collapse=" "))
        if (length(resultTable) == 0) {
            resultTable = currentline
        } else {
            resultTable = rbind(resultTable, currentline)
        }
        li=li+1
    }
    resultTable = data.frame(resultTable)
    return(list(result, resultTable, cnvA2))    
}


#' Reading and sorting scnaFile 
#' 
#'
#' @param scnaFile the SCNA segmentation file
#' @param skipchunk segments with number of data points (probes) no more than this number will be skipped 
#' @return sorted scna segmentation data frame
#' @importFrom dplyr arrange
#' @export
scnaInput <- function(scnaFile, skipchunk=19) {
    #smooth the CN profie according to the minimum segment size (to skip)
    scna = mergeCNA(cnFile = scnaFile, skipchunk = skipchunk)
    #sort by chr and coordinates
    scna = dplyr::arrange(scna, chrom, loc.start, loc.end)
    return(scna)    
}

#' Reading and sorting ssnvFile 
#' 
#'
#' @param ssnvFile the SSNV file
#' @return sorted ssnv data frame
#' @importFrom dplyr arrange
#' @export
ssnvInput <- function(snvFile) {
    if (! file.exists(snvFile)) {
        stop("ssnvFile needs to be provided!")
    }
    ssnv = read.delim(snvFile, stringsAsFactors=FALSE)
    #sort by chr and coordinates
    ssnv = dplyr::arrange(ssnv, chr, pos)
    return(ssnv)
}

#' Generate the SCNA SSNV input for running butte 
#'
#' Four elements will be generated in the output list
#' merged CNA segmentation;
#' SSNV data.frame;
#' cnvHits (index in CNA file, overlapping ssnv)
#' snvHits (index in SSNV file, overlapping with cnv)
#' 
#'
#' @param scnaFile the SCNA segmentation file 
#' @param ssnvFile the SSNV file
#' @param skipchunk segments with number of data points (probes) no more than this number will be skipped 
#' @return list of data input for running butte
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export
cnmutData <- function(scnaFile, ssnvFile, skipchunk = 19) {

    cnvA2 = scnaInput(scnaFile, skipchunk=skipchunk)
    sampAB = ssnvInput(ssnvFile)
    
    #overlap between SCNA and SSNV
    #first make it consistent the chromosome prefix
    cnvSeqNames = as.character(cnvA2$chrom)
    if ( !grepl("chr", cnvA2$chrom[1]) ) {
        cnvSeqNames = paste("chr", cnvA2$chrom, sep="")
    }
    message(cnvSeqNames[1])

    snvSeqNames = sampAB$chr
    if ( !grepl("chr", sampAB$chr[1]) ) {
        snvSeqNames = paste("chr", sampAB$chr, sep="")
    }
    message(snvSeqNames[1])

    #find overlaps
    cnvRangeA = GenomicRanges::GRanges(seqnames = cnvSeqNames,
                                       ranges = IRanges::IRanges(cnvA2$loc.start, end=cnvA2$loc.end),
                                       strand=rep('+',dim(cnvA2)[1]))
    snvRange = GenomicRanges::GRanges(seqnames = snvSeqNames,
                                      ranges = IRanges::IRanges(sampAB$pos, end=sampAB$pos),
                                      strand=rep('+',dim(sampAB)[1]))
    foA = GenomicRanges::findOverlaps(cnvRangeA, snvRange)
    cnvHits = S4Vectors::queryHits(foA)
    snvHits = S4Vectors::subjectHits(foA)
    return(list(cnvA2, sampAB, cnvHits, snvHits))

}



#' Merge the CNA by jumping (neglecting) small segments 
#'
#' When CNA segmentation contains many small segments, one may want to
#' merge the two neighboring segments by skipping the small segment.  
#'
#' @param cnFile the SCNA segmentation file 
#' @param skipchunk segments with number of data points (probes) no more than this number will be skipped 
#' @param correctMale logical, whether or not divide by 2 for the X chromosome (testing)
#' @return data frame of the merged CNA segmentation
#' @export
mergeCNA <- function(cnFile, skipchunk = 19, correctMale = FALSE) {

    if (! file.exists(cnFile)) {
        stop("cnFile (titan segmentation file) not found!")
    }
    message(cnFile)
    cnvA = read.delim(cnFile, stringsAsFactors=FALSE)
    cnvA = cnvA[which(cnvA$num.mark > skipchunk),]                               #skip two few marks
    cp2 <- c(which(cnvA$logcopynumberratio[-1] != cnvA$logcopynumberratio[-nrow(cnvA)] |
                       cnvA$chrom[-1] != cnvA$chrom[-nrow(cnvA)] |
                           cnvA$LOHcall[-1] != cnvA$LOHcall[-nrow(cnvA)]),
             nrow(cnvA))
    cp1 <- c(1,cp2[-length(cp2)]+1)
    cnvA2 <- data.frame(chrom=cnvA$chrom[cp1],
                        loc.start=cnvA$loc.start[cp1],
                        loc.end=cnvA$loc.end[cp2],
                        num.mark=cnvA$num.mark[cp1],         #recal later
                        seg.mean=cnvA$seg.mean[cp1],         #recal later
                        copynumber=cnvA$copynumber[cp1],
                        minor_cn=cnvA$minor_cn[cp1],
                        major_cn=cnvA$major_cn[cp1],
                        allelicratio=cnvA$allelicratio[cp1], #recal later
                        LOHcall=cnvA$LOHcall[cp1],
                        cellularprevalence=cnvA$cellularprevalence[cp1],
                        ploidy=cnvA$ploidy[cp1],
                        normalproportion=cnvA$normalproportion[cp1],
                        logcopynumberratio=cnvA$logcopynumberratio[cp1])
    for (j in 1:length(cp1)) {
        cnvA2$num.mark[j] <- sum(cnvA$num.mark[cp1[j]:cp2[j]])
        cnvA2$seg.mean[j] <- mean(cnvA$seg.mean[cp1[j]:cp2[j]])
        cnvA2$allelicratio[j] <- mean(cnvA$allelicratio[cp1[j]:cp2[j]])
    }
    if (correctMale) {   #testing purpose only, for chromosomeX, divide by 2
        cnvA2[which(as.character(cnvA2$chrom == "X")), "copynumber"] =
            cnvA2[which(as.character(cnvA2$chrom == "X")), "copynumber"]/2
        cnvA2[which(as.character(cnvA2$chrom == "X")), "major_cn"] =
            cnvA2[which(as.character(cnvA2$chrom == "X")), "major_cn"]/2
    }
    return(cnvA2)
}



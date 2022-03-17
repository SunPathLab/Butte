prepareCytoband <- function(cy_file = "/media/ruping/work/stanford_work/OS/timing/cytoband_hg38.txt") {   
    cytoband=read.delim(cy_file, header=F)
    ind = (cytoband[,1] == "chrY" | cytoband[,1] == "chrM")
    cytoband = cytoband[!ind,]
    cytoband[,1] = gsub("chrX","chr23",cytoband[,1])
    cytoband[,1] = as.numeric(sub("chr","",cytoband[,1]))      #how to change chromosome names
    return(cytoband)
}


chrCumCoor <- function(nchr=22, cytoband=NULL,
                       cytofile="/media/ruping/work/stanford_work/OS/timing/cytoband_hg38.txt") {

    if (is.null(cytoband)) {
        cytoband = prepareCytoband(cytofile)
    }
    chrLim=matrix(0,nchr,2)      #chromosome limits
    centromere=numeric()         #centromere limits
    cum=0
    for(i in 1:nchr)
    {
        ind=cytoband[,1]==i
        indP=regexpr("p",cytoband[ind,4])>0
        indQ=regexpr("q",cytoband[ind,4])>0
        
        chrLim[i,]=cbind(min(cytoband[ind,2]),max(cytoband[ind,3]))
        centromere=c(centromere, max(cytoband[ind,3][indP]))
        cum=c(cum,cum[i]+chrLim[i,2])
    }
    return(list(cum,centromere,chrLim))
    
}

cytoband = chrCumCoor(nchr=22, cytofile="../../cytoband_hg38.txt")

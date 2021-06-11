#' @export
patientMutTable <- function(d, pat, norm="_N", titanFolder="./titan/", multiMis=3, lod=4.3, adjv=0.05, correctSampleName=FALSE) {
    samplesall = colnames(d)[grep("maf",colnames(d))]
    samplesall = gsub("maf$", "", samplesall)
    #samplesall = samplesall[grepl("genome", samplesall)]
    
    samples = samplesall[grepl(paste0(pat,"\\D"),samplesall)]
    samplen = samples[grepl(norm,samples)]
    message(paste(samplen, collapse=" "))
    samples = samples[!grepl(norm,samples)]
    message(paste(samples, collapse=" "))
    
    tmp = getSampMutMulti(samples, samplen, d, multiMis, lod)
    tmp = adjust.ccf.titan.multi(tmp, samples, adjv, titanFolder, correctColname=correctSampleName)
    return(tmp)
}

#' @export
write.cnv <- function(samples, titan = "./titan/", outdir="./cnvres/") {
    cnvres <- vector(mode = "list", length = length(samples))
    for(i in 1:length(samples)) {
        cnvresTmp = mergeCNA(titanPath=titan, sn=samples[i])
        if (! dir.exists(outdir)) {
            dir.create(outdir)
        }
        write.table(cnvresTmp, file=paste(outdir, samples[i], ".nclones1.cnv", sep=""),
                    quote=F, row.names=F, sep="\t")
        cnvres[[i]] = cnvresTmp
        names(cnvres)[[i]] = samples[i]
    }
    return(cnvres)
}

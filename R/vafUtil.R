#' Return the allele frequencies of SSNVs for each allele state 
#'
#' Given a SCNA configuration Nt (total copy) and Nb (minor copy),
#' as well as the tumor purity and prevalence of the SCNA,
#' this function returns the expected allele frequencies for possible 
#' allele states.
#'
#' @param ncmut the number of copies for the mutation, default print for "all"
#' @param nt total copy number
#' @param nb copy number of minor allele
#' @param pu tumor purity
#' @param pa prevalence of the SCNA WITHIN the tumor content, default 1 (clonal)
#' @return a named vector of expected allele frequency at each allele state
#' @export
vafEst <- function(ncmut="all", nt, nb, pu, pa=1) {
    vafs = vector()
    epu = pa*pu
    if (ncmut == "all") {
        ncm = 1:(nt-nb)
        for (i in 1:(nt-nb)) {
            vafs = append(vafs, epu*ncm[i]/(epu*nt + 2*(1-epu)))
            names(vafs)[i] = paste0(as.character(ncm[i]),"/",as.character(nt))
        }
    } else {
        vafs = epu*ncmut/(epu*nt + 2*(1-epu))
    }
    return(vafs)
}


#' @export
purityEst <- function(nt, nb, vaf, minor=FALSE, single=FALSE, baf=FALSE) { #suppose major allele
    p = (2*vaf)/((nt-nb)-vaf*(nt-2))
    if (baf) {
        p = (2*vaf-1)/((nt-nb-1)-vaf*(nt-2))
    }
    if (minor) {
        p = (2*vaf)/(nb-vaf*(nt-2))
    }
    if (single) {
        p = (2*vaf)/(1-vaf*(nt-2))
    }
    return(p)
}



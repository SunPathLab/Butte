#' Bootstrap for the CI of timing estimates
#'
#' the corresponding SCNA evolution.
#'
#' @param eventOrdering is the list returned by Butte
#' @param B is the number of bootstrapping samples
#' @param type "bootstrap" for resample the data; or "parametric" to simulate randomized data 
#' @param pi is the estimate of pi you want to use for bootstrapping
#' @param x the data to bootstrap
#' @param m the data to bootstrap
#' @param call list of original argument settings used in calling Butte
#' @return timing estimates of bootstrapped data
#' @export
bootstrapButte <- function(eventOrdering, B, type=c("parametric","bootstrap"),  pi, x, m, call) {

    type<-match.arg(type)

    if (missing(x) | missing(m)) {
        if(!"perLocationProb" %in% names(eventOrdering))
            stop("perLocationProb element is needed if 'x' and 'm' are not supplied.")
        x<-eventOrdering$perLocationProb$x
        m<-eventOrdering$perLocationProb$m
    }
    if(missing(pi)) {
        if(!"pi" %in% names(eventOrdering) | (any(is.na(eventOrdering$pi)) & type=="parametric"))
            stop("invalid values of 'pi' in 'eventOrdering' object")
        pi<-eventOrdering$pi
    }
    if(missing(call)){
        if(!"call" %in% names(eventOrdering))
            stop("invalid values of 'call' in 'eventOrdering' object")
        call<-eventOrdering$call
    }
    
    N <- length(m)
    
    # eventFunctionNames<-c("history","nt","nb","type","qmethod","seqError","purity",
    #                       "minMutations","init","maxiter","tol")
    
    eventFunctionNames<-c("nt","nb","qmethod","seqError","purity",
                          "minMutations","init","maxiter","tol")
    
    requiredNames<-c("alleleSet",eventFunctionNames)
    
    if(!all(requiredNames%in%names(call)))
        stop("Missing call information:",requiredNames[!requiredNames%in%names(call)])
    possAlleles <- call$alleleSet
    A <- call$history[[1]]    #there might be multiple history matrice
    initial <- call$init

    if (type == "parametric") {  #generate random x vector according to m and allele frequencies
        if(ncol(A)!=length(pi))
            stop("'history' from 'call' does not match length of 'pi'")
        q<-as.vector(A%*%pi) 
        q<-q/sum(q)
        if(length(q)!=length(possAlleles))
            stop("'alleleSet' and 'history' from 'call' does not match in size ")
        #deal with any possible problems with numerical inaccuracy
        if(any(q < 0)){
            if(any(q< -1e-10)) stop("Programming error -- negative probabilities")
            else q[q<0]<-0
        } 
        if(any(q > 1)){
            if(any(q> 1+1e-10)) stop("Programming error -- >1 probabilities")
            else q[q>1]<-1
        }
        bootData<-.parametricSimulation(B,m,q,alleleSet=possAlleles,onlyReturnData=TRUE)
    }
    else if (type == "bootstrap") {   #directly sample the data with replacement
        wh <- matrix(sample(1:N, N*B, replace=TRUE), nrow=N, ncol=B)
        bootData <- lapply(1:B, function(kk) {
            cbind(nMutAllele=x[wh[,kk]], nReads=m[wh[,kk]])})
    }
    piBoot<-lapply(bootData, function(z) {
        do.call(Butte, c(list(x=z[,"nMutAllele"],m=z[,"nReads"]), call[eventFunctionNames]))$pi})
    piBoot <- do.call("rbind",piBoot)
    #return(list(piBoot, bootData))
    piBoot[!is.finite(piBoot)] <- NA
    return(piBoot)
}


#assumes alleleSet is the right one (already accounts for sequencing error, normal contamination, etc.)
.parametricSimulation<-function(B, m, q, alleleSet, onlyReturnData=FALSE){
    if(length(q)!=length(alleleSet))
        stop("'alleleSet' and 'q' must be of the same length")
    if(is.null(dim(m))) {   #not a matrix
        nMut <- length(m)
        m <- matrix(m, ncol=B, nrow=nMut, byrow=FALSE)
    }
    else {
        nMut <- nrow(m)
        if(B!=ncol(m)) stop("'B' must be same as number of columns of m")
    }
    if(nMut==0){  #number of mutations equal to zero
        data.frame(AF=rep(NA,nMut),nReads=rep(NA,nMut),nMutAllele=rep(NA,nMut),obsAF=rep(NA,nMut))
    }
    nAlleles<-length(q)
    
    #adjust the allele frequencies to account for Sequencing error:
    
    ###########Random generation of the allele frequencies: #
    #number of each allele frequency will observe
    #a nallele x B matrix of which allele frequency get for each 
    nPerCategory <- rmultinom(n=B, size=nMut, prob=q)
    	
    #repeat each of the alleles according to nPerCategory and then permute them
    #so that each m[i] equally likely to be paired with each allele
    #shuffle (sample without replacement)
    AF<-apply(nPerCategory,2,function(x){sample(rep(alleleSet,times=x))}) #matrix of nMut x B

    #if B=1 of nMut=1, make it into matrix
    if(is.null(dim(AF))){
        AF<-matrix(AF,ncol=B,nrow=nMut,)
    }
    
    #run binomial for each B and turn into reasonable data.frame
    return(lapply(1:B,function(kk){
        x <- rbinom(nMut, size=m[,kk], prob=AF[,kk])
        if(onlyReturnData)
            cbind(nMutAllele=x, nReads=m[,kk])
        else data.frame(AF=AF[,kk],nMutAllele=x,nReads=m[,kk],obsAF=x/m[,kk])
    }))
}

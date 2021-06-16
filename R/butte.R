#' Calculate timing of SCNA
#'
#' This function estimates the timing of SCNAs for the following two scenarios
#' for 2:0 3:1 3:0 4:1, where an identifiable history matrix is available, 
#' this function returns the time estimate and confidence interval for each stage
#' For SCNAs with multiple or non-indentifible matrix, this function instead
#' returns the bounds (lower and upper) of time period for the last stage of
#' the corresponding SCNA evolution.
#'
#' @param x vector of number of reads supporting SSNVs
#' @param m vector of total read depth for SSNVs
#' @param history list of possible histories returned by cnmutHistory
#' @param nt total copy number
#' @param nb copy number of the minor allele
#' @param qmethod suggest to use fullMLE, which is more accurate
#' @param type choose "identifiable" for SCNA 2:0 3:1 3:0 4:1, and "butte" otherwise
#' @param bootstrapCI set to "bootstrap" if non-parametric; or "parametric" 
#' @return A list of possible matrices
#' @export
Butte <- function(x, m, history, nt, nb, qmethod=c("fullMLE","partialMLE"),
                  type=c("identifiable","butte"), seqError=0, bootstrapCI=NULL,
                  B=500, CILevel=0.9, purity=1, verbose=TRUE,
                  returnAssignments=FALSE, minMutations=10, init=NULL,
                  maxiter=100, tol=0.0001, mutationId=1:length(x),...) {

    qmethod <- match.arg(qmethod)
    doCI <- !is.null(bootstrapCI)
    type <- match.arg(type)

    nMuts<-length(x)

    #checking input
    if(length(m)!=nMuts) stop("x and m must be of same length")
    if(length(purity) != 1 | purity > 1 | purity < 0)
        stop("tumor purity must be given and should be single value between 0 and 1.")
    if(length(mutationId)!=length(unique(mutationId)) & verbose)
        warning("mutationId not unique values; using index as id instead.")
    if(length(mutationId)!=nMuts & verbose)
        warning("mutationId an invalid length; using index as id instead.")
    
    if(returnAssignments) returnData <- TRUE
    else returnData <- FALSE
    
    A <- history[[1]]  #use the first history to calculate n steps
    
    K <- ncol(A)-1     #number of steps
    
    #get final allele frequencies
    possAlleles <- vafEst(nt=nt, nb=nb, pu=purity)
    
    #possAlleles <- errorAF(possAlleles,seqError=seqError)

    summaryTable<-cbind(nMuts,length(m))
    colnames(summaryTable)<-c("Total Mutations","Qualify for fitting")
    summaryTable<-t(summaryTable)
    colnames(summaryTable)<-"Number"

    # call arguments for bootstrapping 
    call<-list(alleleSet=possAlleles,
               history=history,
               nt=nt,
               nb=nb,
               type=type,
               qmethod=qmethod,
               seqError=seqError,
               purity=purity,
               minMutations=minMutations,
               init=init,
               maxiter=maxiter, 
               tol=tol)
    if(returnData) inputData=data.frame(mutationId=mutationId,x=x,m=m)
    
    failReason <- NULL   # a sentence tells the reason for failing the calculation
    success <- TRUE      # logical TRUE or FALSE

    # formatting the output list
    dummyOutput<-list("pi"=rep(NA,length=ncol(A)),alleleSet=possAlleles,
                      optimDetails=list(nIter=NA,initial=NA,pathOfQ=NA),
                      perLocationProb=NA,q=NA)
    names(dummyOutput$pi) <- paste("Stage",0:K,sep="")
    # formatting the confidence interval
    ci <- matrix(rep(NA,length=2*length(dummyOutput$pi)),ncol=2)
    row.names(ci) <- names(dummyOutput$pi)
    dummyOutput$piCI <- ci
    
    if ( !is.null(failReason) ) {
        output <- dummyOutput         #failed assign the dummy output
    } else {

        #estimate of time proportion
        piEst <- rep(NA,ncol(A))
        
        # optimizing q vector
        if ( qmethod %in% c("fullMLE") ) {
            ### Estimate the q vector, THE OUTPUT has the SAME elements as dummyOutput
            output<-try(.estimateQ(x, m, possAlleles, history=history[[1]], type=type,
                                   init=init, maxiter=maxiter, tol=tol))
            if( output$optimDetails$nIter == maxiter ) {
                failReason<-paste(failReason,"hit maximum number of iterations",sep=",")
                success<-FALSE
            }
        }
        if( qmethod %in% c("partialMLE") ) { #partialMLE method for estimating q
            mleAlleles <- mleAF(x, m, nt=nt, nb=nb, seqError=seqError, purity = purity)
            output <- try(.estimateQ(x,m, possAlleles, alleleFreq = mleAlleles$alleleSet$Freq,
                                   history=history[[1]], type=type, init=init, maxiter=maxiter, tol=tol))
        }
        if(!inherits(output, "try-error")) {
            # identifiable: directly calculate timing
            if ( type == "identifiable" ) {
                piEst <- solve(A) %*% output$q
                piEst <- as.vector(piEst)
                piEst <- piEst/sum(piEst)
                names(piEst) <- paste("Stage", 0:K, sep="")
                output$pi <- piEst
            }
            
            # non-unique: estimating the bounds instead
            else if ( type == "butte" ) {
                # set to zero for too low mutations
                q2 = output$q
                numq = round(q2 * nMuts)
                q2[which(q2 < 0.01 & numq < 5)] = 0
                q2 = q2/sum(q2)

                buttebounds = try(.lpbounds(q = q2, possible_histories = history))
                piEst <- buttebounds
                piEst[!is.finite(piEst)] <- NA
                names(piEst) <- c("lower", "upper")
                output$pi <- piEst
            }
            
        } else {   #if ERROR occurred when estimating q
            output<-dummyOutput
            failReason<-paste(failReason,"error returned in estimating q",sep=",")
            success<-FALSE
        }
        
        if(!returnAssignments) output<-output[-grep("perLocationProb",names(output))]
        else{ output[["perLocationProb"]]<-data.frame(inputData,output[["perLocationProb"]],
                                                      check.names=FALSE)}
        
        #Calculate Confidence Interval
        if ( doCI ) {
            if ( type == "identifiable" & !any(is.na(output$pi)) ) {   #unique history
                #get pi of bootstrapping samples
                bootsample <- bootstrapButte(call=call, B=B, x=x, m=m,
                                             type=bootstrapCI, pi=output$pi)
                #determining confidence interval
                ci <- t(apply(bootsample, 2, quantile, c((1-CILevel)/2,1-(1-CILevel)/2)))
            }
            else if ( type == "butte" ) {        #non-unique history, get CI for bounds instead
                bootsample <- bootstrapButte(call=call, B=B, x=x, m=m,
                                             type=bootstrapCI, pi=output$pi)
                ci <- t(apply(bootsample, 2, quantile, c((1-CILevel)/2,1-(1-CILevel)/2), na.rm=T))
            }
            else {
                ci<-matrix(rep(NA,length=2*length(output$pi)),ncol=2)
                row.names(ci)<-names(output$pi)
            }
            #colnames(ci)<-c("lower","upper")
            output$piCI <- ci
            #selected elements for output
            output <- output[c("pi","piCI","q","perLocationProb","optimDetails")] 
        } else {
            output<-output[c("pi","q","perLocationProb","optimDetails")]
        }
    }
    return(c(output,list(summaryTable=summaryTable, success=success,
                         failReason=if(!is.null(failReason)) gsub("^,","",failReason)
                                    else failReason,call=call)))
}


#' Use linear programming to solve the bounds of Ta
#'
#' Yunong: add some description here
#'
#' @param q q estimated from data
#' @param possible_histories matrices of possible SCNA-SSNV histories
#' @param scost the cost for slack variables (default 100)
#' @return the lower and upper bounds of the time duration for the last stage
#' @importFrom lpSolve lp
.lpbounds <- function(q, possible_histories, scost=100) {   #modified from Yunong's code with relax  

    lbs = vector()
    ubs = vector()
    nrx = vector()
    for(i in seq(possible_histories)) {

        A = possible_histories[[i]]
        #from A, we can get s vector, which is the copies survive till the end at each time interval
        s = colSums(A)
        #as described by the draft, we can get a matrix M = A-qs', where s' is s transpose, and M*pi = 0
        M = A-q%*%t(s)
        b = numeric(length(q))
        #we should also consider the convexity of the time vector PI
        M = rbind(M, rep(1,ncol(A)))
        b = c(b,1)

        #M add two columns
        relaxCols = diag(dim(M)[1])
        relaxCols = relaxCols[,rep(1:ncol(relaxCols), each=2)]
        relaxCols[,seq(2,dim(relaxCols)[2], by=2)] =  relaxCols[,seq(2,dim(relaxCols)[2], by=2)]*-1
        M = cbind(M, relaxCols)
        
        #notice that, the objective function we want to optimize is
        # f(pi) = [0,0,0,0,0,0 ......,0, 1] * pi, which equals to 
        # the last time interval pi_n
        #define objective function  (represent coefficients in a vector)
        f.obj <- c(rep(0,ncol(A)-1),1)
        f.obj <- c(f.obj, rep(scost, dim(M)[1]*2))  #relax the model
        
        #define constraints
        f.con <- M
        f.dir <- rep("==",nrow(M))
        f.rhs <- b
        
        #message(paste0(i,":\n"))
        lowbound = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
        lb = lowbound$objval - sum(tail(lowbound$solution,dim(M)[1]*2) * scost)
        n_relax = length(which(tail(lowbound$solution,dim(M)[1]*2) > 0))

        f.obj <- c(rep(0,ncol(A)-1),1)
        f.obj <- c(f.obj, rep(-scost, dim(M)[1]*2))  #relax the model
        uppbound = lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs)
        ub = uppbound$objval - sum(tail(uppbound$solution,dim(M)[1]*2) * scost * -1)
        n_relax = n_relax + length(which(tail(uppbound$solution,dim(M)[1]*2) > 0))
        
        if (sum(lowbound$solution) > 0) {
            lbs = append(lbs, lb)
            names(lbs)[length(lbs)] = i
        }
        if (sum(uppbound$solution) > 0) {
            ubs = append(ubs, ub)
            names(ubs)[length(ubs)] = i
            nrx = append(nrx, n_relax)
            names(nrx)[length(nrx)] = i
        }
    }

    #return the min and max only for models with minimum number of non-zero slack variables
    return(c(min(lbs[which(nrx == min(nrx))]), max(ubs[which(nrx == min(nrx))])))
}

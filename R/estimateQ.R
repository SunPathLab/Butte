#' EM algorithm to calculate q
#'
#' Estimate the proportion the proportions of SSNVs at each allele state
#' from sequencing data.
#' 
#'
#' @param x vector of number of reads supporting SSNVs
#' @param m vector of total read depth for SSNVs
#' @param xGreaterZero determines whether the likelihood will account for the fact that only observe X[i]>0
#' @return A list of possible matrices
.estimateQ <- function(x,m,alleleSet,alleleFreq=NULL,history,type=c("identifiable","butte"),
                       init=NULL,maxiter=100, tol=0.0001,xGreaterZero=TRUE, useGradient=FALSE) {
    type = match.arg(type)
        
    N<-length(x)
    possAlleles<-alleleSet
    lengthOfTheta<-length(possAlleles)-1
    if(nrow(history)!=length(possAlleles))
        stop("length of alleleSet must equal the number of rows of history.")

    ##check if it is one of the easy ones:
    Avec<-as.vector(history)
    easyType<-"no"
    if(identical(Avec,c(0, 1, 2, 0))) easyType<-"CNLOH"
    if(identical(Avec,c(1 ,1, 3, 0))) easyType<-"SingleGain"
    #if (easyType != "no") message(easyType)
    
    
    if(is.null(init)) {
        #init<-history%*%rep(1/length(possAlleles),length=length(possAlleles))
	init<-history%*%rep(1/dim(history)[2], length=dim(history)[2])
        init<-as.vector(init)/sum(init)
    }
    
    ##Break it into 2 because useful to have a function that assigns most likely allele to each
    EStep <- function(q) { ####Calculate P(Pi=a | Xi,q) each i and a
        #P(Pi=a | Xi)=P(Xi|Pi=a)P(Pi=a)/sum(P(Xi|Pi=a)P(Pi=a))
        
        #calculate the log of the probability (numerator)
        logPPerAllele<-mapply(possAlleles,q,FUN=function(a,qa){
            dbinom(x=x, size=m, prob=a,log=TRUE)+log(qa) 
        }) #returns a N x numb.Alleles Matrix; needs to be normalized
        
        #-----
        #calculate log of sum of probabilities (denominator)
        #-----
        #note if individual probabilities very small, then denominator sums to zero,
        #so divide by zero -- problem!
        #use log sum exponential trick:
        # log(e[theta1]+e[theta2]+...)=m+log(e[theta1-m]+e[theta2-m]+...), where m=max(theta[i])
        #here theta1=log[P(Xi|Pi=a)]+log[P(Pi=a)]
        #-----

        rowMax<-apply(logPPerAllele,1,max) #max of theta[i]
        logPPerAlleleMinusMax<-sweep(logPPerAllele,1,rowMax,"-") #subtract off m from each row
        logdenom<-rowMax+log(apply(exp(logPPerAlleleMinusMax),1,sum))
        
        #-----
        #calculate log of sum of probabilities (denominator)
        #-----
        #now overall prob=exp[num]/exp[denom]=exp[num-denom]
        #-----
        logPPerAllele<-sweep(logPPerAllele,1,logdenom,"-")
        
	
        # #-----
        # #Old simplistic code:
        # #-----
        #     	PPerAllele<-mapply(possAlleles,q,FUN=function(a,qa){
        # 	dbinom(x=x, size=m, prob=a,log=FALSE)*(qa) 
        # }) #returns a N x numb.Alleles Matrix; needs to be normalized
        # PPerAllele<- prop.table(PPerAllele,1) #divide by the sum of each row
        
        #compare the two:
        # cbind(exp(logPPerAllele),PPerAllele)
        
        return(exp(logPPerAllele))
    }
    
    ###Constraints for the q
    #The feasible region is defined by ui %*% theta - ci >= 0. 
    if(easyType=="no") {
        uiL1<-diag(rep(-1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
        ciL1<-rep(1,lengthOfTheta)
        uiGr0<-diag(rep(1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
        ciGr0<-rep(0,lengthOfTheta)
        if (type == "identifiable") {
            Ainv <- solve(history)
            uiA <- Ainv%*%(rbind(diag(lengthOfTheta),rep(-1,lengthOfTheta)))
            ciA <- - Ainv%*%c(rep(0,lengthOfTheta),1)
            ui <- rbind(uiL1,uiGr0,uiA)
            ci <- c(-ciL1,-ciGr0,ciA)
        } else if (type == "butte") {
            ui <- rbind(uiL1,uiGr0)
            ci <- c(-ciL1,-ciGr0)
        }
    }
    
    #make matrix of (1-a[j])^m[i]
    MStep<-function(PMat=NULL,qinit) {
        #estimate of q constrained so that solve(history)%*%q is all positive
        if(is.null(alleleFreq)) Na<-colSums(PMat)
        #from E-Step, the estimated values from the multinomial...
        #basically just the maximum likelihood estimation if knew Pi exactly
        else Na<-alleleFreq
        if (easyType %in% c("CNLOH","SingleGain")) {
            #print("easytype")
            if (easyType=="CNLOH") { 
                a<-0
                b<-1
            }
            if (easyType=="SingleGain") {
                a<-1/2
                b<-1
            }
            qout<-Na[1]/sum(Na)
            if(qout>b) qout<-b
            if(qout<a) qout<-a
            qout<-c(qout,1-qout)
        } else {
            f<-function(q) {
                qS <- 1-sum(q)   # sum(q) = 1; qS is the last one; constrained
                #-dmultinom(Na, prob=c(q,qS), log = TRUE)
                #just to check; should be equivalent, but faster to not calculated normalizing constant
                ll<- sum(Na*log(c(q,qS)))
                if(xGreaterZero) {
                    #missZero[i]= 1- sum_j{(1-a[j])^m[i] q[j]}
                    missZero<-1-rowSums(mapply(possAlleles,c(q,qS),FUN=function(a,qa){
                        (1-a)^m*qa 
                    }) #returns a N x numb.Alleles Matrix; needs to be normalized
                    )  #sum over all j, so missZero is vector of length i single value for each i
                    ll<-ll-sum(missZero)
                } 
                return( - ll) #because minimizing
            }

            #used in the gradient
            #calculate (1-a[S])^{m[i]} - (1-a[j])^{m[i]} 
            #returns a N x S-1 Matrix; needs to be normalized
            aS <- tail(possAlleles,1)
            # #bug reported
            # num<-mapply(head(possAlleles,-1),q,FUN=function(a,qa){
            # 	(1-aS)^m - (1-a)^m
            # })
            #fix with
            num<-sapply(head(possAlleles,-1),FUN=function(a,qa) {
                (1-aS)^m - (1-a)^m
            })
            gr<-function(q) {
                qS<-1-sum(q)
                NS<-tail(Na,1)
                g<-(head(Na,-1)/q - NS/qS)
                #was - (head(Na,-1)/q - NS/qS*q) #which is error; why got wrong answer.
                if(xGreaterZero) {
                    #missZero[i]= 1- sum_j{(1-a[j])^m[i] q[j]}
                    missZero<-1-rowSums(mapply(possAlleles,c(q,qS),FUN=function(a,qa){
                        (1-a)^m*qa 
                    }) #returns a N x numb.Alleles Matrix; needs to be normalized
                    )  #sum over all j, so missZero is vector of length i single value for each i
                    g<-g-colSums(sweep(num,1,missZero,"/"))
                }
                return(-g) #because minimizing
            }
            theta <- head(qinit,-1)
            if ( lengthOfTheta==1 ) {
                upper<-min((ci/ui)[ui<0])
                lower<-max((ci/ui)[ui>0])
                out<-optimize(f=f,interval=c(lower,upper),tol=.Machine$double.eps)
                qout<-c(out$minimum,1-out$minimum)
            } else {
                out<-constrOptim(theta=theta,f = f, grad=if(useGradient) gr else NULL, ui=ui, ci=ci)
                #when use gr, get wrong answer -- returns me to init, doesn't iterate, etc. WHY????
                qout<-c(out$par,1-sum(out$par))      # return the q estimates
            }
        }
        names(qout)<-names(possAlleles)
        return(qout)
    }
    
    #print("gr uses - (head(Na,-1)/q - NS/qS) ")
    TotalStep<-function(q) {
        MStep(EStep(q),qinit=q)
    }
    if(is.null(alleleFreq)) {
        qOld <- init;
        qNew<-TotalStep(qOld)
        allQ<-cbind(qOld,qNew)
        nIter<-1
        #convMessages<-c(firstOut$convergence)
        while(sum(abs(qNew - qOld)) >= tol & nIter <= maxiter) {
            qOld<-qNew
            pMat<-EStep(qOld)
            qNew<-MStep(pMat,qinit=qOld)
            #qNew<-c(mstepOut$par,1-sum(mstepOut$par))
            #convMessages<-c(convMessages,mstepOut$convergence)
            nIter <- nIter + 1
            allQ<-cbind(allQ,qNew)
            #if(nIter %% 20 ==0) print(sprintf("finished iteration %s", nIter))
        }
        pMat<-EStep(qNew)
    }
    #given the alleleFreq obtained from partialMLE, optimize q directly
    else {
        qNew<-MStep(qinit=init)
        pMat<-NA
        nIter<-NA
        allQ<-NA
    }
    names(qNew)<-names(possAlleles)
    #return optimization results and details
    return(list(q=qNew,perLocationProb=pMat,alleleSet=possAlleles,
                optimDetails=list(nIter=nIter,initial=init,pathOfQ=allQ)))
}


mleAF <- function(x, m, nt, nb, seqError=0, purity=1) {

    if(length(x)!=length(m)) stop("x and m must be the same length")

    maxCopy = nt-nb
    r <- 1:maxCopy 

    alleles <- vafEst(nt=nt, nb=nb, pu=purity)
    alleles <- errorAF(alleles, seqError=seqError)
    
    out <- sapply(alleles, function(z){dbinom(x, size=m, prob=z)})
    colnames(out) <- r
    assign <- r[apply(out, 1, which.max)]
    assignAlleles <- factor(alleles[apply(out,1,which.max)], levels=alleles)
    afdf <- data.frame(tumorAF=r/nt, AF=alleles, table(assignAlleles))
    row.names(afdf) <- NULL
    afdf <- afdf[,-3]
    
    outList<-list(perLocationProb=out,assignments=data.frame(nCopies=assign,nt=nt,tumorAF=assign/nt))
    outList<-c(outList,list(alleleSet=afdf))
    return(outList)
    
}


errorAF <- function(trueAF, seqError=0) {
    ((3-4*seqError)*trueAF+seqError)/(3-2*seqError)
}

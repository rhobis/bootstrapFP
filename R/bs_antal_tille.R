### Antal and Till? bootstrap for pps sampling ---------------------------------

#' Antal and Till? (2011) bootstrap for unequal probability sampling
#' 
#' Draw B bootstrap resamples according to Antal and Till? (2011) direct
#' bootstap method for unequal probability sampling.
#' Note that this method do not need a double bootstrap.
#'
#'@param ys values of the variable of interest for the original sample
#'@param pks vector of first-order inclusion probabilities for sampled units
#'@param B integer scalar, number of bootstrap resamples to draw from the pseudo-population
#'@param smplFUN function used to select the sample according to the
#'original sampling design. The function must accept the vector of inclusion
#'probabilities as input and return a vector o 0's and 1's (as the functions of
#'the library sampling).
#'@param approx_method method used to approximate Dkk with function \code{Dkk_aprox()}
#'@param ... added to ignore useless arguments
#'
#'@return a list of two elements, a vector of K average bootstrap totals and
#'a vector of K variance estimates.
#' 
#' 
#' @references 
#' 
#' Antal, E.; Till?, Y., 2011. A Direct Bootstrap Method for Complex Sampling
#' Designs From a Finite Population. Journal of the American Statistical Association, 106:494, 534-543,
#' doi: 10.1198/jasa.2011.tm09767
#' 
#' Antal, E.; Till?, Y., 2014. A new resampling method for sampling designs without
#' replacement: the doubled half bootstrap. Computational Statistics, 29(5), 1345-1363.
#' doi: 10.10007/s00180-014-0495-0



bs_antal_tille <- function(ys, pks, B, smplFUN, approx_method = c("Hajek", "DevilleTille"), ...) {
    
    ### Initial values ---
    approx_method <- match.arg(approx_method)
    n    <- length(ys)
    Dkk  <- Dkk_approx(pks, n, method=approx_method)
    phi  <- (1 - Dkk)
    n1   <- sum(phi)
    case <- n - n1
    
    ht <- vector('numeric', length=B)
    
    if(case >= 2){
        ### CASE 1 (algorithm 4) -----------------------------------------------
        
        
        n1_int <- as.integer(n1)
        if( n1 != n1_int ){
            q <- n1_int - n1 + 1
            phi1 <- decompose_phi(phi)
            rand_ind <- sample( x=1:2, size=1, prob=c(q, 1-q) )
            m <- phi1[[ rand_ind ]][[2]]
            phi1 <- phi1[[ rand_ind ]][[1]]
        }
        
        ### Bootstrap ---
        for( b in seq_len(B) ){
            #step 1
            S <- smplFUN(phi1)
            ind_not_S <- which( S == 0 )
            #step 2
            oos <- one_one( n = length(ind_not_S) ) 
            #step 3
            S[ind_not_S] <- S[ind_not_S] + oos
            
            # compute HT estimate
            ht[b] <- sum(S * ys/pks)
        }#end for
        
    } else {
        ### CASE 2 (algorithm 5) -----------------------------------------------   
        
        #step 1
        psi <- ( 1 - sampling::inclusionprobabilities(1-phi, n = 2) ) 
        q   <- case/2
        
        for( b in seq_len(B)){
            u   <- runif(1)
            if( u <= q){ #step 2
                S <- smplFUN(psi)
                oos <- one_one( n = 2 )
                ind_not_S <- which( S == 0)
                S[ ind_not_S ] <- S[ ind_not_S ] + oos
            } else { #step 3
                S <- rep(1, length=n) 
            }
            ht[b] <- sum(S * ys/pks)
        }
    }
    
    ### Return results ---
    return( list( Tb = mean(ht),
                  Vb = var(ht) )
    )
    
}


### ----------------------------------------------------------------------------
#' Approximation of the Dkk quantity 
#' 
#' see Antal and Till? (2011)
#' 
#' @param pks vector of first-order inclusion probabilities for sampled units
#' @param n integer, the sample size
#' @param method a string indicating which approximation method is to be used
#' 
#' @return a numeric vector with approximated Dkk

Dkk_approx <- function(pks, n, method=c("Hajek", "DevilleTille")) {
    
    method <- match.arg(method)
    
    if(method == "Hajek"){
        ck <- (n/(n-1)) * (1-pks) 
        Dkk <- ck - (ck*ck / sum(ck))
    } else if(method == "DevilleTille" ){
        if(n>2){
            Dkk <- (1-pks)    
        } else stop("Deville and Tille's method does not converge for n=2")
    } 
    
    return(Dkk)
}


### ----------------------------------------------------------------------------
#' Decompose phi vector
#' 
#' Decompose vector phi in a convex combination of two vectors phi1 and phi2, 
#' such that the sum of $phi1_i$ is the integer part of phi and the sum of $phi2_i$
#' is the integer part of phi plus 1 (see Antal and Till? (2011) bootstrap procedure
#' for unequal probability sampling, p. 539 - Algorithm 4, Case 1) 
#' 
#' @param phi vector of inclusion probabilities for Antal and Till? (2011) bootstrap, 
#' given by 1 - D_kk
#' 
#' @return a list with the two vectors in which \code{phi} is decomposed

decompose_phi <- function(phi) {
    ### Initial values ---
    nphi   <- sum(phi)
    m1     <- as.integer(nphi)
    m2     <- m1 + 1
    u      <- runif(length(phi)) #random vector with sum != 0
    
    ### Decompose vector phi
    phi1   <- phi + (m1 - nphi) * u/sum(u)
    phi2   <- phi + (m2 - nphi) * u/sum(u)
    
    ### Return results ---
    return( list(phi1 = list(phi = phi1, m = m1),
                 phi2 = list(phi = phi2, m = m2))
    )
}


### ----------------------------------------------------------------------------
#' Select a one-one sampling 
#' 
#' A one-one sampling is a design for which the random variables Sk, representing
#' the number of times unit k is present in the sample, have expectation and
#' variance equal to 1. Proposed by Antal and Till? (2011, 2014).
#' 
#' @param n integer, the sample size
#' @param method algorithm to be used, currently, only the doubled half sampling
#' is implemented (Antal and Till?, 2014). See Details section.
#' 
#' @details Antal and Till? proposed two procedures that lead to one-one samplings.
#' The first one (Antal and Till?, 2011a) in more complex and makes use of a simple 
#' random Sampling with over-replacement (Antal and Till?, 2011b), 
#' and it is called by setting \code{method = "over-replacement"};
#' the second one (Antal and Till?, 2014) is the doubled half sampling, which is 
#' simpler and quickier to compute, it can be called by setting
#' \code{method = "doubled-half"} and is the default option.
#' 
#' @return an integer vector of size \code{n}, indicating how many times each unit is
#' present in the sample
#' 
#' 

one_one <- function(n, method = c("doubled-half", "over-replacement") ){
    method <- match.arg(method)
    
    if(method == "doubled-half"){
        S <- doubled_half( n )
    }else if(method == "over-replacement"){
        if( n == 2 ){
            S <- vector( 'integer', length = n )
            ind <- sample( 1:2, size=1)
            S[ind] <- 2
        }else{
            #step 1
            m <- as.integer( 0.5 * ( 1 + sqrt ( (4*n^2 + 5*n -1) / (n-1) ) ) )
            alpha <- ( m*(n-1)*(m+1) - n*(n+1) ) / ( 2*m*(n-1) )
            #step 2
            nt <- sample( c(m, m+1), size=1, prob = c(alpha, 1-alpha))
            #step 3
            sB <- over_replacement( N=n, n=nt) #vector of length N=n
            #step 4
            sA <- drop( rmultinom(1, size = (n-nt), prob = rep(1/n, n)) )
            #step 5
            S <- sA + sB
        }
    }else stop("Incorrect choice of method.")
    
    return( S ) 
}


### ----------------------------------------------------------------------------
#' Select a doubled-half sampling (Antal and Till?, 2014)
#' 
#' @param n integer scalar representing sample size
#' 
#' @return an integer vector of size \code{n}, indicating how many times each unit is
#' present in the sample


doubled_half <- function( n ){
    if(n<2) stop("sample size must be > 1")
    S <- vector( 'integer', length=n )
    
    ### Select sample ---
    if (n %% 2 == 0){
        ind <- sample(1:n, size = n/2)       
        S[ind] <- 2
    }else {
        ind <- sample(1:n, size = (n-1)/2)       
        S[ind] <- 2
        u <- runif(1)
        if( u <= 0.25){
            ind2 <- sample( ind, size = 1)
        }else{
            ind2 <- sample( which(S==0), size = 1)
        }
        S[ind2] <- S[ind2]+1
    }
    ### Return selected sample
    return( S )
}


### ----------------------------------------------------------------------------
#' Select a simple random sampling with over-replacement 
#' 
#' Used for resampling procedures. Proposed by Antal and Till? (2011).
#' 
#' @param N integer, the population size
#' @param n integer, the sample size
#' 
#' @return an integer vector of size n, indicating how many times each unit is
#' present in the sample
#' 
#' @references 
#' Antal, E.; & Till?, Y. (2011). Simple random sampling with over-replacement. 
#' Journal of Statistical Planning and Inference, 141(1), 597-601.

over_replacement <- function( N, n ){
    if (N<2 | n<2) stop("Population and sample size must be > 1")
    
    ### Initialize variables ---
    sumS <- 0
    nk <- n
    S <- vector('integer', length=N)
    
    ### Sample selection ---
    for( k in seq_len(N)){
        j <- 0:nk
        p <- choose(N-k-1+nk-j, nk-j) / choose(N-k+nk, nk)
        Sk <- sample(j, size=1, prob=p)
        sumS <- sumS + Sk
        nk <- n - sumS
        S[k] <- Sk
    }
    
    ### Return selected sample ---
    return( S )
}

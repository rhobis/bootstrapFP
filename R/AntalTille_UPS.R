#' Antal and Tillé (2011) Bootstrap for Unequal Probability Sampling without replacement
#' 
#' Draw B bootstrap samples according to Antal and Tillé (2011) direct
#' bootstap method for Unequal Probability Sampling.
#' Note that this method does not need a double bootstrap.
#'
#' @param ys values of the variable of interest for the original sample
#' @param pks vector of first-order inclusion probabilities for sampled units
#' @param B integer scalar, number of bootstrap resamples to draw from the pseudo-population
#' @param smplFUN a function that takes as input a vector of length N of 
#' inclusion probabilities and return a vector of length N, either logical or a 
#' vector of 0s and 1s,  where \code{TRUE} or \code{1} indicate sampled
#' units and \code{FALSE} or \code{0} indicate non-sample units.
#' @param approx_method method used to approximate the variance Dkk.
#'
#' @return a list of two elements, a vector of K average bootstrap totals and
#' a vector of K variance estimates.
#' 
#' 
#' 
#' 
#' @references 
#' 
#' Antal, E.; Tillé, Y., 2011. A Direct Bootstrap Method for Complex Sampling
#' Designs From a Finite Population. Journal of the American Statistical Association, 106:494, 534-543,
#' doi: 10.1198/jasa.2011.tm09767
#' 
#' Antal, E.; Tillé, Y., 2014. A new resampling method for sampling designs without
#' replacement: the doubled half bootstrap. Computational Statistics, 29(5), 1345-1363.
#' doi: 10.10007/s00180-014-0495-0
#' 
#' 
#'
#' @keywords internal



AntalTille2011_ups <- function(ys, pks, B, smplFUN, approx_method = c("Hajek", "DevilleTille")) {
    
    ### Initialisation ---
    approx_method <- match.arg(approx_method)
    n    <- length(ys)
    
    # Approximation of the Dkk quantity, see Antal and Tille' (2011)
    if(identical(approx_method, "Hajek") ){
        ck  <- (n/(n-1)) * (1-pks) 
        Dkk <- ck - ck^2/sum(ck)
    } else if( identical(approx_method, "DevilleTille") ){
        #(a solution does not always exist)
        Dkk <- 1-pks
    } 
    
    phi  <- 1-Dkk
    n1   <- sum(phi)
    case <- n-n1   #this value determines the algorithm to use (<2 or >=2)
    
    ht <- vector('numeric', length=B)
    
    if(case >= 2){
        ### CASE 1 (algorithm 4) -----------------------------------------------
        
        phi <- define_phi(phi)
        
        ### Bootstrap ---
        for( b in seq_len(B) ){
            #step 1
            S <- as.logical(smplFUN(phi))
            #step 2
            ind_not_S <- which( S == 0 )
            oos <- one_one( n = length(ind_not_S) ) # one-one sampling
            #step 3
            S[ind_not_S] <- oos
            
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
                S <- as.logical(smplFUN(psi))
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
    return( var(ht) )
    
}




### ----------------------------------------------------------------------------
#' Define the phi vector
#' 
#' Define the phi vector used to select the first sample in Antal & Tillé (2011)
#' bootstrap (algorithm 4, first step).
#' If the sum of the elements of \eqn{\phi}{phi} is not an integer, phi is decomposed 
#' in a convex combination of two vectors \eqn{\phi_1}{phi1} and \eqn{\phi_2}{phi2}, 
#' such that the sum of \eqn{\phi_1i}{phi1_i} is the integer part of \eqn{\sum phi_i}{sum(phi_i)} 
#' and the sum of \eqn{\phi_2i}{phi2_i} is the integer part of \eqn{\sum phi_i}{sum(phi_i)} plus 1 
#' [see Antal and Tille' (2011) bootstrap procedure
#' for unequal probability sampling, p. 539 - Algorithm 4, Case 1]
#' The procedure used to decompose the vector \eqn{\phi}{phi} is described in the 
#' answer to this question: https://math.stackexchange.com/questions/2700483/vector-decomposition-into-a-convex-combination-of-two-vectors-with-constraints-o
#' 
#' @param phi vector of inclusion probabilities for Antal and Tillé (2011) bootstrap, 
#' given by 1 - D_kk
#' 
#' @return a list with the two vectors in which \code{phi} is decomposed
#' 
#' 
#' @keywords internal

define_phi <- function(phi) {
    ### Initial values ---
    np <- sum(phi)
    
    if( is_wholenumber(np)) {
        
        return(phi)
        
    } else {
        
        m1 <- as.integer(np)
        m2 <- m1 + 1
        # u <- runif(length(phi)) #random vector with sum != 0
        u  <- rep(1, length(phi))
        
        ### Decompose vector phi
        ur    <- u/sum(u)
        phi1  <- phi + (m1 - np)*ur
        phi2  <- phi + (m2 - np)*ur
        
        ### Randomly draw either phi1 or phi2
        q    <- m1 - np + 1
        r    <- sample( x=1:2, size=1, prob=c(q, 1-q) )
        out  <- if(r==1) phi1 else phi2
        
        return( out )
        
    }
}


### ----------------------------------------------------------------------------
#' Select a one-one sampling 
#' 
#' A one-one sampling is a design for which the random variables Sk, representing
#' the number of times unit k is included in the sample, have expectation and
#' variance equal to 1. Proposed by Antal and Tille' (2011, 2014).
#' 
#' @param n integer, the sample size
#' @param method algorithm to be used, either doubled half sampling or srs with over-replacement. 
#' See the Details section.
#' 
#' @details Antal and Tillé proposed two procedures that lead to one-one samplings.
#' The first one (Antal and Tillé, 2011a) in more complex and makes use of a simple 
#' random Sampling with over-replacement (Antal and Tillé, 2011b), 
#' and it is called by setting \code{method = "over-replacement"}.
#' The second one (Antal and Tillé, 2014) is the doubled half sampling, which is 
#' simpler and quickier to compute, and can employed by setting
#' \code{method = "doubled-half"}; this is the default option.
#' 
#' @return an integer vector of size \code{n}, indicating how many times each unit is
#' present in the sample
#' 
#' 
#' 
#' @importFrom stats rmultinom
#' 
#' 
#' @keywords internal

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
#' Select a doubled-half sampling (Antal and Tille', 2014)
#' 
#' @param n integer scalar representing sample size
#' 
#' @return an integer vector of size \code{n}, indicating how many times each unit is
#' present in the sample
#' 
#' @keywords internal


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
#' Used for resampling procedures. Proposed by Antal and Tille' (2011).
#' 
#' @param N integer, the population size
#' @param n integer, the sample size
#' 
#' @return an integer vector of size n, indicating how many times each unit is
#' present in the sample
#' 
#' @references 
#' Antal, E.; Tillé, Y. (2011). Simple random sampling with over-replacement. 
#' Journal of Statistical Planning and Inference, 141(1), 597-601.
#' 
#' 
#' @keywords internal

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

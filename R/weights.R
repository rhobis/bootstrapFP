#' Bootstrap with Adjusted Weights 
#'
#' Compute bootstrap estimates according to Bootstrap Weights procedures by
#' Rao et Al. (1992) and Chipperfield and Preston (2007).
#'
#'@param ys values of the variable of interest for the original sample
#'@param N scalar, representing the population size
#'@param B integer scalar, number of bootstrap resamples to draw from the pseudo-population
#'@param method a string indicating the bootstrap method to be used; 
#'       available methods are "RaoWuYue" and "ChipperfieldPreston".
#'@param ... added to ignore useless arguments
#'
#'@return a list of two elements, a vector of K average bootstrap totals and
#'a vector of K variance estimates.
#'
#'
#'@references
#'
#'Rao J. N. K.; Wu C. F. J.; Yue K. (1992).  Some recent work on resampling methods 
#'for complex surveys. Journal of the American Statistical Association, 83(401), 620-630.
#'
#'Chipperfield J.; Preston J. (2007).Efficient bootstrap for business surveys.
#'Survey Methodology, 33(2), 167-172.
#'
#'
#'
#' @keywords internal





bootstrap_weights <- function(ys, 
                              N, 
                              B, 
                              method = c("RaoWuYue", "ChipperfieldPreston")
){
    
    method <- match.arg(method)
    
    ### Initialize variables ---
    n  <- length(ys)
    f  <- n/N
    w  <- N/n #original weights
    
    if( identical(method, "RaoWuYue")){
        nn <- n-1
        smplFUN <- function() sampling::srswr(n=nn, N=N)
        den <- n-1
    }else if( identical(method, "ChipperfieldPreston")){
        nn <- floor(n/2)
        smplFUN <- function() sampling::srswor(n=nn, N=N)
        den <- n - nn
    }
    
    tot <- vector('numeric', length=B)
    
    ### Bootstrap procedure ---
    for(b in seq_len(B)){
        
        #Generate weights
        m  <- smplFUN()
        a  <- 1 + sqrt( nn*(1-f)/den ) * (n*m/nn - 1)
        
        ws <- a*w
        
        #Estimate total
        tot[b] <- sum(ws*ys)
    }
    
    
    ### Return results ---
    Tot <- sum(w*ys) # estimator in the original sample
    return( sum( (tot-Tot)^2 )/B ) #Variance estimation as in Rao et Al (1992)  
    
}





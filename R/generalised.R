#' Generalised Bootstrap 
#'
#' Compute bootstrap estimates according to Generalised Bootstrap procedure by
#' Beaumont and Patak (2012)
#'
#'@param ys values of the variable of interest for the original sample
#'@param pks inclusion probabilities for units in the sample
#'@param B integer scalar, number of bootstrap resamples to draw from the pseudo-population
#'@param distribution the distribution from which to generate the weights 
#'adjustments. One of \code{uniform}, \code{normal} or \code{lognormal}.
#'
#'@return a list of two elements, a vector of K average bootstrap totals and
#'a vector of K variance estimates.
#'
#'
#'@references
#'
#'Bertail, P., & Combris, P. (1997). Bootstrap généralisé d'un sondage. 
#'Annales d'Economie et de Statistique, 49-83.
#'
#'Beaumont, J. F., & Patak, Z. (2012). On the generalized bootstrap for sample 
#'surveys with special attention to Poisson sampling. 
#'International Statistical Review, 80(1), 127-148.
#'
#'
#' @importFrom stats rnorm runif
#' 
#' 
#' @keywords internal





generalised <- function(ys, 
                           pks, 
                           B, 
                           distribution = c("uniform", "normal", "exponential", "lognormal") 
                        ){
    
    distribution <- match.arg(distribution)
    
    ### Initialize variables ---
    n  <- length(ys)
    ht <- vector('numeric', length=B)
    w  <- 1/pks #original weights    
    
    ### Bootstrap procedure ---
    for(b in seq_len(B)){
        ### Generate weights
        switch(distribution,
               "uniform"     = {  q <- sqrt( 3*(1-pks) )
                                  a <- runif(n, 1-q, 1+q)
                               },
               "normal"      = {  a <- rnorm(n, 1, 1-pks) },
               "exponential" = {  a <- 1 + ( exp(n)-1 ) * sqrt(1-pks)},
               "lognormal"   = {  q <- log(2-pks)
                                  a <- exp( rnorm(n, -0.5*q, q) )
                               }
        )
        ws   <- a*w
        ht[b] <- sum(ys*ws)
    }
    
    HT <- drop( sampling::HTestimator(ys, pks) ) 
    ### Return results ---
    # return( list( Tb = mean(ht), 
    #               Vb = (sum((ht-HT))^2)/B ) #Variance estimation as in Beaumont and Patak(2012)
    return( (sum((ht-HT))^2)/B ) #Variance estimation as in Beaumont and Patak(2012) ) 
    
}



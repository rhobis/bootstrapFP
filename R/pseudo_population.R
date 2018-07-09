#' Pseudo-population bootstrap for simple random sampling
#'
#' @param y vector of sample values
#' @param N scalar, represents the population size
#' @param B scalar, number of bootstrap replications
#' @param D scalar, number of replications for the double bootstrap (when applicable)
#' @param method a string indicating the bootstrap method to be used, available
#'     methods are: 'Gross', 'Booth', 'ChaoLo85', 'ChaoLo94', 'BickelFreedman', 'Sitter'
#'
#'

ppBS_srs <- function(y, N, B, D=1, method) {
    
    ### Check input ---
    method <- match.arg(method,
                        c('Gross',
                          'Booth',
                          'ChaoLo85',
                          'ChaoLo94',
                          'BickelFreedman',
                          'Sitter')
    )
    
    if( identical(method, 'Gross') & !is_wholenumber(N/n)) 
        stop("Gross method can be used only when N/n is an integer. Please, choose another method!")
    
    
    ### Initialise quantities ---
    n <- length(y)
        
    # quantities for the random population, Uc
    if( identical(method, 'Gross') ){
        k <- floor( N/n ) 
        select_Uc <- function( ... ){
            return( NULL )
        }
        
    } else if( identical(method, 'Booth') ){
        
        k  <- floor( N/n )
        select_Uc <- function(){
            s <- as.logical( sampling::srswor( (N-n*k), n) )
            return( y[s] )
        }
        
    } else if( identical(method, 'ChaoLo94') ){
        
        k <- floor( N/n )
        select_Uc <- function(){
            s <- as.logical( sampling::srswr( (N-n*k), n) )
            return( y[s] )
        }
        
    } else if( identical(method, 'ChaoLo85') ){
        
        k <- floor( N/n )
        
        G <- function(x){
            return( (1 - n/x) * x*(n-1) / ((x-1)*n) )
        }
        q <- (G(N) - G(n*(k+1)) ) / (G(n*k) - G( n*(k+1)) )
        select_Uc <- function(){
            if( q <= runif(1) ){
                return( NULL )
            } else{
                return( y )
            } 
        }
        
    } else if( identical(method, 'BickelFreedman') ){
        
        k <- floor( N/n )
        q <- ( 1 - (N-n*k)/n ) * ( 1 - (N-n*k)/(N-1))
        select_Uc <- function(){
            if( q <= runif(1) ){
                return( NULL )
            } else{
                return(y)
            } 
        }
        
    } else if( identical(method, 'Sitter') ){
        
        f  <- n/N
        
        a1 <- ( n*k - n + 1 ) / ( n*(n-1)*(n*k - 1) )
        a2 <- k / ( n * ( n*(k+1) - 1 ) )
        q  <- (1-f)/(n*(n-1)) - a2
        q  <- q / ( a1-a2 )
        
        select_Uc <- function(){
            if( q <= runif(1) ){
                return( NULL )
            } else{
                return( y )
            } 
        }
        
    } 
        
        
    
    Uf  <- rep(ys, each = k )

    Vb <- vector('numeric', length = D)
    for(d in seq_len(D) ){

        
        Uc <- select_Uc()
        n1 <- ifelse(is.null(Uc) & identical(method, 'Sitter'), n-1, n)
        U  <- c( Uf, Uc)
        
        tot <- vector('numeric', length = B)
        ### Bootstrap procedure
        for(b in seq_len(B)){
            #select resample
            sb <- U[ sample(length(U), n1, replace = FALSE) ]
            tot[b] <- N * mean(sb)
        }
        Vb[d] <- var(tot)
    }
    
    ### Return results ---
    return( mean(Vb) )
    
}



ppBS_ups <- function(y, pik, B, D) {
    
    
    
    
}
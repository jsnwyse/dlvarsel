###### Default Prior for Distributed Lag parameters

dlprior <- function( p, q, m, l, weight.scheme ) {
  
  prior <- list()
  
  sigsq <- 1
  prior$precision_static <- diag( 1/c(4, rep( sigsq, q-1 )) )
  prior$precision_dynamic <- diag( 1/sigsq, m * p )
  
  if( weight.scheme == "b-spline" )
  {
    prior$prior.mean.theta <- rep(NULL,m)
    prior$prior.sd.theta <- rep(NULL,m)
    prior$prior.rate.theta <- rep(NULL,m)
  }
  if( weight.scheme == "exp-b-spline" )
  {
    prior$prior.mean.theta <- rep( rep(0,m), p ) 
    prior$prior.sd.theta <- rep( rep(1,m), p )
    prior$prior.rate.theta <- rep(NULL,2*p)
  }
  if( weight.scheme == "nealmon" ) 
  {
    prior$prior.mean.theta <- rep( l/10, 2*p ) 
    prior$prior.sd.theta <- rep( (l/10)/4 , 2*p )
    prior$prior.rate.theta <- rep(10,2*p)
  }
  
  return(prior)
}

dlstart <- function( p, m, prior, weight.scheme ) {
  
  if ( p == 0 ) { return(NULL) }

  if( weight.scheme == "b-spline" ) return(NULL)
  
  if( weight.scheme == "exp-b-spline" ) init.theta <- rnorm( p * m, 0, 1 )

  if( weight.scheme == "nealmon" ) 
  {
    init.theta <- rnorm( 2*p, prior$prior.mean.theta, prior$prior.sd.theta )
    idx <- seq(2, 2*p, by=2)
    init.theta[ idx ] <- - rexp( length(idx), rate = prior$prior.rate.theta[idx] )
  }
  
  return( list(init.theta=init.theta) )  
}


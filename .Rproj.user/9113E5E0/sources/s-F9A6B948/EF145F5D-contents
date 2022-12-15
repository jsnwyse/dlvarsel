dlvarsel <- function( y, W = NULL, X = NULL, formula = NULL, tau, l, 
                      weight.scheme = c("exp-b-spline","b-spline","nealmon"), m=NULL , 
                      prior = NULL, start = NULL,
                      nsamp = 10^4, burn = 10^3, thin = 1, adapt.interval = 50, sampmeth = c("lang","smmala"), dl.varsel = FALSE,
                      verbose = FALSE )
{
  
  if( is.null(formula) & !is.null(W) ) formula <- as.formula(paste( "~ ", paste( names(W), collapse = "+")))
  W <- if( is.null(W) ) matrix( rep(1,length(y)), nrow=length(y), ncol=1) else model.matrix( formula, data = W )  # gives intercept by default
  n <- length(y)
  q <- max(1,ncol(W))
  p <- max(1,nrow(X)/l)
  if( is.null(weight.scheme) ) weight.scheme <- "exp-b-spline"
  if( is.null(m) )
  {
    m <- if( weight.scheme == "nealmon" ) 1 else 4
  }else{
    if( m < 1 ) stop("m number of basis functions must be 2 or greater")
    if( weight.scheme == "nealmon" & m > 1 )
    {
      warning("m must equal 1 for nealmon weighting- this has been reset automatically")
      m <- 1 
    }
  }
  static.vars <- ifelse( is.null(W), FALSE, TRUE )
  if( is.null(X) ) stop("not possible to run without distributed lag variables")
  dl.vars <- TRUE 

  if( length(y) != n ) stop("length of y is not the same as number of rows in W")
  
  if( is.null(prior) ) prior <- dlprior( p, q, m, l, weight.scheme )
  if( is.null(start) ) start <- dlstart( p, m, prior, weight.scheme )
  
  basis.eval <- if( weight.scheme == "nealmon" ) cbind( (1:l), (1:l)^2 ) else bs( 1:l, degree=m )
  ord.nealmon <- if( weight.scheme == "nealmon" ) 2 else m
  if( weight.scheme == "exp-b-spline" ) m <- 1

  gammasamp <- numeric( nsamp * q )
  betasamp <- numeric( nsamp * m * p )
  etasamp <- numeric( nsamp * n )
  thetasamp <- numeric( nsamp * p * ord.nealmon )
  varsamp <- numeric( nsamp * (q + p) )
  loglik <- numeric( nsamp )
  prprobsamp <- numeric( nsamp )
  lvn.stepsize <- c(0.75,0.0)
  acc.rates <- numeric(2)

  if( is.null(sampmeth) ) sampmeth <- "lang"
  langevin.samp <- 2 - match( sampmeth, c("lang","smmala") )

  if( weight.scheme == "nealmon" ) ws <- 0
  if( weight.scheme == "exp-b-spline" ) ws <- 1
  if( weight.scheme == "b-spline" ) ws <- 2
  weight.scheme <- ws
  
  if( dl.varsel & ( weight.scheme == 0 || weight.scheme == 2) ) stop("variable selection option only available for 'exp-b-spline' weighting scheme")
  
  var.ind <- if( dl.varsel ) c( rep(1,q), rep(0,p)) else rep(1, q+p)
  sd.varsel <- 2.5
  var.ind.prior.prob <- rep(0.5, q+p)
  hprior.var.ind <- c(1,1)

  out <- .C( "distlag_MCMC_bspl", as.integer(n), as.integer(y),
             as.integer(q), as.double(W), 
             as.integer(p), as.double(X),
             as.integer(l), as.integer(m), 
             as.double(tau), as.double(prior$precision_static),
             as.double(prior$precision_dynamic), as.double( rep(basis.eval,p) ) ,
             gammasamp = as.double(gammasamp), etasamp = as.double(etasamp),
             betasamp = as.double(betasamp), thetasamp = as.double(thetasamp),
             varsamp = as.integer(varsamp), 
             loglik = as.double(loglik), prprobsamp = as.double(prprobsamp),
             as.integer(nsamp), as.integer(burn),
             as.integer(thin), as.integer(static.vars),
             as.integer(dl.vars), as.integer(weight.scheme),
             as.double(start$init.theta), 
             as.double(prior$prior.mean.theta), as.double(prior$prior.sd.theta), as.double(prior$prior.rate.theta),
             as.integer(ord.nealmon),
             sd.theta.prop = as.double(lvn.stepsize), as.integer(adapt.interval),
             accrt = as.double(acc.rates), as.integer(langevin.samp), as.integer(dl.varsel),
             as.integer(var.ind), as.double(sd.varsel), as.double(var.ind.prior.prob),
             as.double( hprior.var.ind ), as.integer(verbose),
             PACKAGE="dlvarsel" )
  
  ret <- list(samples=list())
  ret$samples$gamma <- matrix( out$gammasamp, nrow=q )
  ret$samples$beta <- matrix( out$betasamp, nrow=m*p )
  ret$samples$theta <- matrix( out$thetasamp, nrow=ord.nealmon*p )
  ret$samples$eta <- matrix( out$etasamp, nrow=n )
  ret$samples$varsamp <- matrix( out$varsamp, nrow= q+p )
  ret$samples$loglik <- out$loglik
  ret$samples$priorprob <- out$prprobsamp
  ret$langstep <- out$sd.theta.prop[1]
  ret$acc.rt <- out$accrt
  names( ret$acc.rt ) <- c("dl pars", "var sel")
  
  return(ret)
  
}
GATE.Gibbs <- function( TestStatistic, burnin, thin, size, alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE ){

    G <- length( TestStatistic )
    mgs <- array(0, G)
    pi1.post <- array(0, size)
    pi2.post <- array(0, size)
    mu.post <- array(0, size)
    sigmaSq.post <- array(0, size)

    for(g in 1:G)
        mgs[g] <- TestStatistic[[g]]$mg
    
    ## Initialization
    pi1.tmp <- rbeta(1, alpha1, beta1)
    pi2.tmp <- rbeta(1, alpha2, beta2)
    mu.tmp <- 2
    sigmaSq.tmp <- 1
    
    for( num.rep in 1:( burnin + thin*size ) )
    ## for(num.rep in 2:100)
    {
        alpha1_prime <- alpha1; beta1_prime <- beta1; alpha2_prime <- alpha2; beta2_prime <- beta2;
            
        ## Step one, update theta.g, update theta.jg
        ## ########################################
        mu.post.mean <- 0 ## This is the variable for calculating the posterior distribution of mu, \sum_{gj} \theta_{gj}X_{gj}
        total.no.zero <- 0 ## This is the quantity:  \sum_{gj}\theta_{gj}
        sigma.post.scale <- 0 ## This is the quantity: \sum_{gj}(x_{gj}-mu)^2 
        
        ## ## Update theta.g, and theta.jg
        for(g in 1:G)
        {
            ## ## Update theta.g
            TestStatistic[[g]]$f0x = dnorm( TestStatistic[[g]]$X ) * 10
            TestStatistic[[g]]$f1x = dnorm( TestStatistic[[g]]$X, mu.tmp, sqrt(sigmaSq.tmp) )*10
            TestStatistic[[g]]$f0prod <- prod( TestStatistic[[g]]$f0x )
            
            post.prob <- (1-pi1.tmp) * TestStatistic[[g]]$f0prod/( (1-pi1.tmp)*TestStatistic[[g]]$f0prod + pi1.tmp* ( prod( (1-pi2.tmp)*TestStatistic[[g]]$f0x + pi2.tmp*TestStatistic[[g]]$f1x) - (1-pi2.tmp)^(mgs[g]) * TestStatistic[[g]]$f0prod )/( 1- (1-pi2.tmp)^(mgs[g]) ) )
            TestStatistic[[g]]$theta.g.tmp <- ( runif(1) > post.prob )          
                  
            ## ## Update theta.jg
            TestStatistic[[g]]$theta.gj.tmp = array(0, mgs[g])
            if( TestStatistic[[g]]$theta.g.tmp == 1){
              Ratio <- 1 - ( (1-pi2.tmp)*TestStatistic[[g]]$f0x )/( (1-pi2.tmp)*TestStatistic[[g]]$f0x+pi2.tmp*TestStatistic[[g]]$f1x )
              tmp.tmp <- runif( mgs[g] )
              if( sum(tmp.tmp < Ratio) > 0 )
                {
                  TestStatistic[[g]]$theta.gj.tmp <- ( tmp.tmp < Ratio )
                }else{
                  TestStatistic[[g]]$theta.gj.tmp <- array(0, mgs[g] )
                  TestStatistic[[g]]$theta.gj.tmp[ which(tmp.tmp == min(tmp.tmp))] <- 1
              }  
            }

            mu.post.mean <- mu.post.mean + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * TestStatistic[[g]]$X ) ## This is the variable for calculating the posterior distribution of mu, \sum_{gj} \theta_{gj}X_{gj}
            total.no.zero <- total.no.zero + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp ) ## This is the quantity:  \sum_{gj}\theta_{gj}
            sigma.post.scale <- sigma.post.scale + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * (TestStatistic[[g]]$X-mu.tmp)^2) ## This is the quantity: \sum_{gj}(x_{gj}-mu)^2 
  
            
            alpha1_prime <- alpha1_prime + TestStatistic[[g]]$theta.g.tmp
            beta1_prime <- beta1_prime + 1 - TestStatistic[[g]]$theta.g.tmp
            alpha2_prime <- alpha2_prime + TestStatistic[[g]]$theta.g.tmp* sum( TestStatistic[[g]]$theta.gj.tmp)
            beta2_prime <- beta2_prime + TestStatistic[[g]]$theta.g.tmp * (mgs[g]-sum(TestStatistic[[g]]$theta.gj.tmp) )
          
        }
        ## Update pi1.tmp and pi2.tmp
        pi1.tmp <- rbeta(1, alpha1_prime, beta1_prime)
        pi2.tmp <- rbeta(1, alpha2_prime, beta2_prime)

        ## Update mu.tmp and sigmaSq.tmp
        if( total.no.zero > 0)
        {
            shrinkage.factor <- (1/prec.mu)/( 1/prec.mu + sigmaSq.tmp/total.no.zero)
            mu.tmp <- rnorm(1, shrinkage.factor * mu.post.mean/total.no.zero, shrinkage.factor * sigmaSq.tmp/total.no.zero )
            sigmaSq.tmp <- 1/( rgamma(1, shape= shape.sigma + total.no.zero/2, scale=1/( 1/scale.sigma + sigma.post.scale/2 ) ) )
        }
        if( total.no.zero == 0)
        {
            mu.tmp <- 2
            sigmaSq.tmp <- 1
        }

        
        if( ((num.rep-burnin)%%thin ==0 ) && (num.rep>burnin) )
        {
            kk = (num.rep-burnin)/thin
            pi1.post[kk] <- pi1.tmp
            pi2.post[kk] <- pi2.tmp
            mu.post[kk] <- mu.tmp
            sigmaSq.post[kk] <- sigmaSq.tmp
        }
        if( verbose==TRUE)
        if( num.rep%%1000==0)
            print( paste("Total number of rep: ", burnin+thin*size, "; Current rep: ", num.rep, sep="" ) )

      }


    
    TestStatistic <- GATE.localfdr.Bayes(TestStatistic, median(pi1.post), median(pi2.post), median(mu.post), median(sigmaSq.post) )
    list( TestStatistic=TestStatistic, mu.post=mu.post, pi1.post=pi1.post, pi2.post=pi2.post, sigmaSq.post=sigmaSq.post)
    
}


GATE.localfdr.Bayes <- function(TestStatistic, pi1, pi2.1, mu, sigmaSq)
{
  pi0 <- 1-pi1
  pi2.0 <- 1-pi2.1
  G <- length( TestStatistic )
  for( g in 1:G )
    {
        mg <- TestStatistic[[g]]$mg
        TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X )
        TestStatistic[[g]]$f1x <- dnorm( TestStatistic[[g]]$X, mu, sqrt(sigmaSq) )
        TestStatistic[[g]]$fx <- pi2.0 * TestStatistic[[g]]$f0x + pi2.1 * TestStatistic[[g]]$f1x ## Calcualte the marginal density for each x_{gj}
        ## Temparary quantities
        f0.prod <- prod(  10 * TestStatistic[[g]]$f0x )
        f.prod <- prod(  10 * TestStatistic[[g]]$fx )
        f1.prod <- prod( 10 * TestStatistic[[g]]$f1x )
        ## Calcualte the group-local fdr
        TestStatistic[[g]]$fdr.g <- pi0 * f0.prod/( pi0*f0.prod + pi1 * ( f.prod - pi2.0^mg * f0.prod )/( 1- pi2.0^mg) )
        ## Calculate the within group local fdr fdr_{j|g}
        den <- f.prod - pi2.0^mg * f0.prod
        num <- pi2.0 * TestStatistic[[g]]$f0x/TestStatistic[[g]]$fx * f.prod - pi2.0^mg * f0.prod
        TestStatistic[[g]]$fdr.j.g <- num/den
    }
  TestStatistic
}


GATE.emGibbs <- function( TestStatistic, burnin, thin, size, alpha1=1, beta1=1, alpha2=1, beta2=1, verbose=TRUE ){

    G <- length( TestStatistic )
    mgs <- array(0, G)
    pi1.post <- array(0, size)
    pi2.post <- array(0, size)
    
    ## Initialization
    pi1.tmp <- rbeta(1, alpha1, beta1)
    pi2.tmp <- rbeta(1, alpha2, beta2)
    min.X <- 0
    max.X <- 0
        
    for( g in 1:G)
    {
        mgs[g] <- TestStatistic[[g]]$mg
        TestStatistic[[g]]$theta.g.tmp <- ( runif(1) < pi1.tmp )
        TestStatistic[[g]]$theta.jg.tmp <- ( runif(mgs[g] ) < pi2.tmp )
        TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X ) * 10 ## Times 10 to avoid to small numbers
        TestStatistic[[g]]$f0prod <- prod( TestStatistic[[g]]$f0x )
        TestStatistic[[g]]$f1x.post <- array( 0, mgs[g] )
        
        if( min.X > min(TestStatistic[[g]]$X) )
            min.X <- min( TestStatistic[[g]]$X )
        if( max.X < max( TestStatistic[[g]]$X) )
            max.X <- max( TestStatistic[[g]]$X)
        
    }

    for( num.rep in 1:( burnin + thin*size ) )
    ##    for(num.rep in 1:100)
    {
        ## Step one, estimate the kernel denstity, update theta.g, update theta.jg
        ## ########################################
        ## ## Calculate kernel density
        x_all_kernel <- c()
        for(g in 1:G)
        {
            if( TestStatistic[[g]]$theta.g.tmp== 1)
                x_all_kernel <- c(x_all_kernel, TestStatistic[[g]]$X[ which( TestStatistic[[g]]$theta.jg.tmp==1)] )
         }
        x_all_kernel <- c(x_all_kernel, min.X, max.X )
        
        ker.den <- approxfun( density( x_all_kernel ) )

        alpha1_prime <- alpha1; beta1_prime <- beta1; alpha2_prime <- alpha2; beta2_prime <- beta2;
        ## ## Update theta.g, and theta.jg
        for(g in 1:G)
        {
            ## ## Update theta.g
            TestStatistic[[g]]$f1x = ker.den( TestStatistic[[g]]$X )*10
            
            post.prob <- (1-pi1.tmp) * TestStatistic[[g]]$f0prod/( (1-pi1.tmp)*TestStatistic[[g]]$f0prod + pi1.tmp* ( prod( (1-pi2.tmp)*TestStatistic[[g]]$f0x + pi2.tmp*TestStatistic[[g]]$f1x) - (1-pi2.tmp)^(mgs[g]) * TestStatistic[[g]]$f0prod)/( 1- (1-pi2.tmp)^(mgs[g]) ) )
            TestStatistic[[g]]$theta.g.tmp <- ( runif(1) > post.prob )

          
            ## ## ## Update theta.jg
            ## TestStatistic[[g]]$theta.jg.tmp = array(0, mgs[g])
            ## if( TestStatistic[[g]]$theta.g.tmp == 1){
            ##     Ratio <- ( TestStatistic[[g]]$f1x/TestStatistic[[g]]$f0x * pi2.tmp/( 1-pi2.tmp) )
            ##     if( max(Ratio) >= 1){
            ##         TestStatistic[[g]]$theta.jg.tmp[ which( Ratio >= 1) ] <- 1
            ##     }else{
            ##         TestStatistic[[g]]$theta.jg.tmp[ which(Ratio==max(Ratio))] <- 1
            ##     }
            ## }
            
            ## ## Update theta.jg
            TestStatistic[[g]]$theta.jg.tmp = array(0, mgs[g])
            if( TestStatistic[[g]]$theta.g.tmp == 1){
              Ratio <- 1 - ( (1-pi2.tmp)*TestStatistic[[g]]$f0x )/( (1-pi2.tmp)*TestStatistic[[g]]$f0x+pi2.tmp*TestStatistic[[g]]$f1x )
              tmp.tmp <- runif( mgs[g] )
              if( sum(tmp.tmp < Ratio) > 0 )
                {
                  TestStatistic[[g]]$theta.gj.tmp <- ( tmp.tmp < Ratio )
                }else{
                  TestStatistic[[g]]$theta.gj.tmp <- array(0, mgs[g] )
                  TestStatistic[[g]]$theta.gj.tmp[ which(tmp.tmp == min(tmp.tmp))] <- 1
              }  
            }

            
            alpha1_prime <- alpha1_prime + TestStatistic[[g]]$theta.g.tmp
            beta1_prime <- beta1_prime + 1 - TestStatistic[[g]]$theta.g.tmp
            alpha2_prime <- alpha2_prime + TestStatistic[[g]]$theta.g.tmp* sum( TestStatistic[[g]]$theta.gj.tmp)
            beta2_prime <- beta2_prime + TestStatistic[[g]]$theta.g.tmp * (mgs[g]-sum(TestStatistic[[g]]$theta.gj.tmp) )
        }
        ## Update pi1.tmp and pi2.tmp
        pi1.tmp <- rbeta(1, alpha1_prime, beta1_prime)
        pi2.tmp <- rbeta(1, alpha2_prime, beta2_prime)

        if( ((num.rep-burnin)%%thin ==0 ) && (num.rep>burnin) )
        {
            kk = (num.rep-burnin)/thin
            pi1.post[kk] <- pi1.tmp
            pi2.post[kk] <- pi2.tmp
            for(g in 1:G)
            {
                TestStatistic[[g]]$f1x.post <- TestStatistic[[g]]$f1x.post + TestStatistic[[g]]$f1x
            }
        }
        if( verbose==TRUE)
        if( num.rep%%1000==0)
            print( paste("Total number of rep: ", burnin+thin*size, "; Current rep: ", num.rep, sep="" ) )

    }

    
    for( g in 1:G)
        TestStatistic[[g]]$f1x.post <- TestStatistic[[g]]$f1x.post/size

    TestStatistic <- GATE.localfdr.nonpar(TestStatistic, mean(pi1.post), mean(pi2.post) )

    list( TestStatistic=TestStatistic, pi1.post=pi1.post, pi2.post=pi2.post)
    
}

GATE.localfdr.nonpar <- function(TestStatistic, pi1, pi2.1)
{
  pi0 <- 1-pi1
  pi2.0 <- 1-pi2.1
  G <- length( TestStatistic )
  for( g in 1:G )
    {
        mg <- TestStatistic[[g]]$mg
        TestStatistic[[g]]$f0x = dnorm( TestStatistic[[g]]$X ) * 10
        
        TestStatistic[[g]]$fx <- pi2.0 * TestStatistic[[g]]$f0x + pi2.1 * TestStatistic[[g]]$f1x.post ## Calcualte the marginal density for each x_{gj}
      ## Temparary quantities
        f0.prod <- prod(  TestStatistic[[g]]$f0x )
        f.prod <- prod(  TestStatistic[[g]]$fx )
        f1.prod <- prod( TestStatistic[[g]]$f1x )
        ## Calcualte the group-local fdr
        TestStatistic[[g]]$fdr.g <- pi0 * f0.prod/( pi0*f0.prod + pi1 * ( f.prod - pi2.0^mg * f0.prod )/( 1- pi2.0^mg) )
        ## Calculate the within group local fdr fdr_{j|g}
        den <- f.prod - pi2.0^mg * f0.prod
        num <- pi2.0 * TestStatistic[[g]]$f0x/TestStatistic[[g]]$fx * f.prod - pi2.0^mg * f0.prod
        TestStatistic[[g]]$fdr.j.g <- num/den
    }
  TestStatistic
}


GATE.localfdr <- function(TestStatistic, pi1, pi2.1, L, muL, sigmaL, cL)
{
  pi0 <- 1-pi1
  pi2.0 <- 1-pi2.1
  G <- length( TestStatistic )
  for( g in 1:G )
    {
      mg <- TestStatistic[[g]]$mg
      TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X ) ## Calculate the density under the null for each x_{gj}
      TestStatistic[[g]]$f1x <- unlist( lapply( TestStatistic[[g]]$X, function(x, L, muL, sigmaL, cL){ sum( cL * dnorm(x, muL, sigmaL) ) }, L, muL, sigmaL, cL ) ) ## Calculate the density function udner the alternative for each x_{gj}
      TestStatistic[[g]]$fx <- pi2.0 * TestStatistic[[g]]$f0x + pi2.1 * TestStatistic[[g]]$f1x ## Calcualte the marginal density for each x_{gj}
      ## Temparary quantities
      f0.prod <- prod( 10* TestStatistic[[g]]$f0x )
      f.prod <- prod( 10* TestStatistic[[g]]$fx )
      f1.prod <- prod( 10* TestStatistic[[g]]$f1x )
      ## Calcualte the group-local fdr
      TestStatistic[[g]]$fdr.g <- pi0 * f0.prod/( pi0*f0.prod + pi1 * ( f.prod - pi2.0^mg * f0.prod )/( 1- pi2.0^mg) )
      ## Calculate the within group local fdr fdr_{j|g}
      den <- f.prod - pi2.0^mg * f0.prod
      num <- pi2.0 * TestStatistic[[g]]$f0x/TestStatistic[[g]]$fx * f.prod - pi2.0^mg * f0.prod
      TestStatistic[[g]]$fdr.j.g <- num/den
      ## Calculate the posterior probability taht both theta and theta_{j|g}=1
      TestStatistic[[g]]$prob.theta.1.theta.j.g.1 <- (1 - TestStatistic[[g]]$fdr.j.g)*( 1-TestStatistic[[g]]$fdr.g ) 

      ## Calculate the posterior prob that theta=1, theta_{j|g}=1, and x_{jg} is generated from the l-th component.
      TestStatistic[[g]]$m.j.g <- array(0, c(mg, L) )
      for( j in 1:mg )
        {
          TestStatistic[[g]]$m.j.g[j, ] <- cL * dnorm( TestStatistic[[g]]$X[j], muL, sigmaL )
          TestStatistic[[g]]$m.j.g[j, ] <- TestStatistic[[g]]$m.j.g[j,] / sum( TestStatistic[[g]]$m.j.g[j, ] )
          TestStatistic[[g]]$m.j.g[j, ] <- TestStatistic[[g]]$m.j.g[j,] * TestStatistic[[g]]$prob.theta.1.theta.j.g.1[j]
        }
    }
  TestStatistic
}



GATE.em <- function( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN=FALSE)
## em <- function( TestStatistic, L=2, pi1.ini=0.5, pi2.1.ini=0.5, muL.ini=rnorm(L), sigmaL.ini=rgamma(L, 2, 2), cL.ini=rdirichlet(1, array(1,L), DELTA=0.001 )
{
  G <- length( TestStatistic )
  mgs <- array(0, G)
  pi1.old <- pi1.ini;  pi2.1.old <- pi2.1.ini;  muL.old <- muL.ini;  sigmaL.old <- sigmaL.ini;  cL.old <- cL.ini ## Set the initial value as the old iteration
  pi1.new <- pi1.old;  pi2.1.new <- pi2.1.old;  muL.new <- muL.old;  sigmaL.new <- sigmaL.old;  cL.new <- cL.old ## Set the new iteration as the old iteration
  delta <- 1 ## Set the

  while( delta> DELTA){
    ##while( itr < 7 ){
    TestStatistic <- GATE.localfdr(TestStatistic, pi1.old, pi2.1.old, L, muL.old, sigmaL.old, cL.old)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
      total.sigmaL <- array(0, L )
      Part.A <- 0
      Part.B <- 0
      fdr.g.all <- array(0, G)
      
      for(g in 1: G)
      {
          fdr.g.all[g] <- TestStatistic[[g]]$fdr.g
          Part.A <- Part.A + (1-fdr.g.all[g]) * sum( TestStatistic[[g]]$fdr.j.g )
          Part.B <- Part.B + (1-fdr.g.all[g]) * sum( 1-TestStatistic[[g]]$fdr.j.g )
          
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
      ## Update all the parameters
      pi1.new <- 1 - total.fdr.g/G
      ## To update pi2.1.new, we use the Newton Raphlson Algorithm
      pi2.1.new <- solve.pi2.1.new(Part.A, Part.B, fdr.g.all, mgs )
           
      
      cL.new <- total.m.j.g.L/total.theta.1.theta.j.g.1
      muL.new <- total.muL/total.m.j.g.L
      if( sigma.KNOWN==FALSE ){
          for(g in 1:G)
              total.sigmaL <- total.sigmaL + apply( ( TestStatistic[[g]]$X %x% array(1, c(1,L)) - array(muL.new, c(1,L))%x% array(1, c(TestStatistic[[g]]$mg, 1) ))^2 * TestStatistic[[g]]$m.j.g, 2, sum )
          sigmaL.new <- sqrt( total.sigmaL/total.m.j.g.L )
      }else{
          sigmaL.new <- sigmaL.old
      }
      ## Calculate the squared distance of all the parameters
      ## delta <- max( abs(pi1.new-pi1.old), abs(pi2.1.new-pi2.1.old), abs(muL.new-muL.old), abs(sigmaL.new-sigmaL.old), abs(cL.new-cL.old ))
      delta <- (pi1.new-pi1.old)^2 + (pi2.1.new-pi2.1.old)^2 + sum( (muL.new-muL.old)^2 ) + sum( (cL.new-cL.old)^2 ) + sum( (sigmaL.new-sigmaL.old)^2 )
      ## Update the iteration by setting the current iteration as the old iteration
      pi1.old <- pi1.new; pi2.1.old <- pi2.1.new; muL.old <- muL.new; sigmaL.old <- sigmaL.new; cL.old <- cL.new;
      print(delta)
  }

  em.esti <- list( pi1=pi1.new, pi2.1=pi2.1.new, muL=muL.new, sigmaL=sigmaL.new, cL=cL.new, L )
  em.esti
}

solve.pi2.1 <- function(Part.A, Part.B, fdr.g.all, mgs )
{
    x.old <- 0.5
    x.new <- 0.7
    G <- length( fdr.g.all )
    DELTA.NR <- 0.0001
    while( abs(x.old-x.new) > DELTA.NR )
    {
        x.old <- x.new
        f.x0 <- -Part.A/(1-x.old) + Part.B/x.old - sum( mgs^2*(1-fdr.g.all)* (1-x.old)^(mgs-1)/(1-(1-x.old)^(mgs +(mgs==1) -1) )  + (mgs==1)*( 1-fdr.g.all )/x.old )
        ##        f.prime.x0 <- -Part.A/(x.old-1)^2 - Part.B/x.old^2 + sum( mgs^2*(1-fdr.g.all)*( 1-x.old)^(mgs-2)*(mgs-1+(1-x.old)^mgs)/( 1-(1-x.old)^mgs)^2 )
        f.prime.x0 <- -Part.A/(x.old-1)^2 - Part.B/x.old^2 + sum( (mgs^2*(1-fdr.g.all)*( 1-x.old)^(mgs-2)*(mgs-1+(1-x.old)^mgs)/( 1-(1-x.old)^mgs)^2 )*(mgs!=1) + (mgs==1)*(1-fdr.g.all)/x.old^2 )
        x.new <- x.old - f.x0/f.prime.x0
        if( x.new > 0.95)
            x.new <- 0.95
        if( x.new < 0.05 )
            x.new <- 0.05
        print(x.new)
    }
    x.new
}

solve.pi2.1.new <- function(Part.A, Part.B, fdr.g.all, mgs )
{
    ## This replaces  the Newton Raphlson Algorithm for solving pi2.1 by searching the solution
    x0 <- c(1:900)/1001+0.1
    y0 <- x0
    for( i in 1:900)
    {
        y0[i] <- Part.A * log( 1-x0[i] ) + Part.B * log( x0[i] ) - sum( log( 1-(1-x0[i])^mgs) * mgs * (1-fdr.g.all ) )
    }
    x0[ which(y0==max(y0) ) ]
}



GATE.em.old <- function( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN=FALSE)
## em <- function( TestStatistic, L=2, pi1.ini=0.5, pi2.1.ini=0.5, muL.ini=rnorm(L), sigmaL.ini=rgamma(L, 2, 2), cL.ini=rdirichlet(1, array(1,L), DELTA=0.001 )
{
  G <- length( TestStatistic )
  mgs <- array(0, G)
  pi1.old <- pi1.ini;  pi2.1.old <- pi2.1.ini;  muL.old <- muL.ini;  sigmaL.old <- sigmaL.ini;  cL.old <- cL.ini ## Set the initial value as the old iteration
  pi1.new <- pi1.old;  pi2.1.new <- pi2.1.old;  muL.new <- muL.old;  sigmaL.new <- sigmaL.old;  cL.new <- cL.old ## Set the new iteration as the old iteration
  delta <- 1 ## Set the

  while( delta> DELTA){
    ##while( itr < 7 ){
    TestStatistic <- GATE.localfdr(TestStatistic, pi1.old, pi2.1.old, L, muL.old, sigmaL.old, cL.old)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
      total.sigmaL <- array(0, L )
      Part.A <- 0
      Part.B <- 0
      fdr.g.all <- array(0, G)
      
      for(g in 1: G)
      {
          fdr.g.all[g] <- TestStatistic[[g]]$fdr.g
          Part.A <- Part.A + (1-fdr.g.all[g]) * sum( TestStatistic[[g]]$fdr.j.g )
          Part.B <- Part.B + (1-fdr.g.all[g]) * sum( 1-TestStatistic[[g]]$fdr.j.g )
          
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
      ## Update all the parameters
      pi1.new <- 1 - total.fdr.g/G
      pi2.1.new <- Part.B/(Part.A+Part.B)
           
      
      cL.new <- total.m.j.g.L/total.theta.1.theta.j.g.1
      muL.new <- total.muL/total.m.j.g.L
      if( sigma.KNOWN==FALSE ){
          for(g in 1:G)
              total.sigmaL <- total.sigmaL + apply( ( TestStatistic[[g]]$X %x% array(1, c(1,L)) - array(muL.new, c(1,L))%x% array(1, c(TestStatistic[[g]]$mg, 1) ))^2 * TestStatistic[[g]]$m.j.g, 2, sum )
          sigmaL.new <- sqrt( total.sigmaL/total.m.j.g.L )
      }else{
          sigmaL.new <- sigmaL.old
      }
      ## Calculate the squared distance of all the parameters
      delta <- max( abs(pi1.new-pi1.old), abs(pi2.1.new-pi2.1.old), abs(muL.new-muL.old), abs(sigmaL.new-sigmaL.old), abs(cL.new-cL.old ))
      ## Update the iteration by setting the current iteration as the old iteration
      pi1.old <- pi1.new; pi2.1.old <- pi2.1.new; muL.old <- muL.new; sigmaL.old <- sigmaL.new; cL.old <- cL.new;
      ## print(delta)
  }

  em.esti <- list( pi1=pi1.new, pi2.1=pi2.1.new, muL=muL.new, sigmaL=sigmaL.new, cL=cL.new, L )
  em.esti
}



GATE.em2 <- function( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN=FALSE)
## em <- function( TestStatistic, L=2, pi1.ini=0.5, pi2.1.ini=0.5, muL.ini=rnorm(L), sigmaL.ini=rgamma(L, 2, 2), cL.ini=rdirichlet(1, array(1,L), DELTA=0.001 )
{
  G <- length( TestStatistic )
  mgs <- array(0, G)
  pi1.old <- pi1.ini;  pi2.1.old <- pi2.1.ini;  muL.old <- muL.ini;  sigmaL.old <- sigmaL.ini;  cL.old <- cL.ini ## Set the initial value as the old iteration
  pi1.new <- pi1.old;  pi2.1.new <- pi2.1.old;  muL.new <- muL.old;  sigmaL.new <- sigmaL.old;  cL.new <- cL.old ## Set the new iteration as the old iteration
  delta <- 1 ## Set the

  while( delta> DELTA){
    ##while( itr < 7 ){
    TestStatistic <- GATE.localfdr(TestStatistic, pi1.old, pi2.1.old, L, muL.old, sigmaL.old, cL.old)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
    total.sigmaL <- array(0, L )
    for(g in 1: G)
      {
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
    ## Update pi2.1.new and cL.new
    pi2.1.new <- total.theta.1.theta.j.g.1 / total.theta.g.1
    cL.new <- total.m.j.g.L/total.theta.1.theta.j.g.1


    TestStatistic <- GATE.localfdr(TestStatistic, pi1.old, pi2.1.new, L, muL.old, sigmaL.old, cL.new)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
    total.sigmaL <- array(0, L )
    for(g in 1: G)
      {
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
    ## Update pi1
    pi1.new <- 1 - total.fdr.g/G


    TestStatistic <- GATE.localfdr(TestStatistic, pi1.new, pi2.1.new, L, muL.old, sigmaL.old, cL.new)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
    total.sigmaL <- array(0, L )
    for(g in 1: G)
      {
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
    ## Update muL
    muL.new <- total.muL/total.m.j.g.L

    
    if( sigma.KNOWN==FALSE ){
      TestStatistic <- GATE.localfdr(TestStatistic, pi1.new, pi2.1.new, L, muL.new, sigmaL.old, cL.new)
      total.fdr.g <- 0
      total.theta.g.1 <- 0
      total.theta.1.theta.j.g.1 <- 0
      total.m.j.g.L <- array(0, L)
      total.muL <- array(0, L)
      total.sigmaL <- array(0, L )
      for(g in 1: G)
        {
          mgs[g] <- TestStatistic[[g]]$mg
          total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
          total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
          total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
          total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
          total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
        }
      for(g in 1:G)
        total.sigmaL <- total.sigmaL + apply( ( TestStatistic[[g]]$X %x% array(1, c(1,L)) - array(muL.new, c(1,L))%x% array(1, c(TestStatistic[[g]]$mg, 1) ))^2 * TestStatistic[[g]]$m.j.g, 2, sum )
      sigmaL.new <- sqrt( total.sigmaL/total.m.j.g.L )
    }else{
      sigmaL.new <- sigmaL.old
    }
    ## Calculate the squared distance of all the parameters
    delta <- max( abs(pi1.new-pi1.old), abs(pi2.1.new-pi2.1.old), abs(muL.new-muL.old), abs(sigmaL.new-sigmaL.old), abs(cL.new-cL.old ))
    ## Update the iteration by setting the current iteration as the old iteration
    pi1.old <- pi1.new; pi2.1.old <- pi2.1.new; muL.old <- muL.new; sigmaL.old <- sigmaL.new; cL.old <- cL.new;
    print(delta)
  }

  em.esti <- list( pi1=pi1.new, pi2.1=pi2.1.new, muL=muL.new, sigmaL=sigmaL.new, cL=cL.new, L )
  em.esti
}


GATE.1 <- function(TestStatistic, alpha=0.05 )
  {
    G <- length( TestStatistic )
    all.lfdr <- c()
    for(g in 1:G)
      {
        TestStatistic[[g]]$lfdr.temp <- 1 - ( 1-TestStatistic[[g]]$fdr.g )*( 1-TestStatistic[[g]]$fdr.j.g)
        all.lfdr <- c( all.lfdr, TestStatistic[[g]]$lfdr.temp )
      }
    all.lfdr <- sort( all.lfdr, decreasing=FALSE )
    total.Rej <- max ( ( cumsum( all.lfdr ) < c(1:length(all.lfdr))*alpha ) * c(1:length(all.lfdr) ) )

    if( total.Rej ==0 ){
      for(g in 1:G)
        {
          TestStatistic[[g]]$gate1.bg.rej <- 0
          TestStatistic[[g]]$gate1.wg.rej <- ( TestStatistic[[g]]$fdr.j.g)*0
          TestStatistic[[g]]$lfdr.temp <- NULL
        }
    }else{
      thresh <- all.lfdr[ total.Rej ]
      for(g in 1:G)
        {
          TestStatistic[[g]]$gate1.wg.rej <- ( TestStatistic[[g]]$lfdr.temp <= thresh )
          TestStatistic[[g]]$gate1.bg.rej <- ( sum(TestStatistic[[g]]$gate1.wg.rej) > 0 )
          TestStatistic[[g]]$lfdr.temp <- NULL
        }
    }
    TestStatistic
  }


GATE.2 <- function(TestStatistic, alpha=0.05, eta=alpha)
  {
    G <- length(TestStatistic)
    group.fdr.star.g <- array(0, G)
    Rgs <- array(0, G)
    ## Go in to each group, calculate the potential number of rejections and mark the corresponding hypothesis, calculate group.fdr.star.g
    for(g in 1:G)
      {
        mg <- TestStatistic[[g]]$mg
        order.fdr.j.g <- sort( TestStatistic[[g]]$fdr.j.g, decreasing=FALSE )
        R_g <- max ( ( cumsum( order.fdr.j.g)/c(1:mg) <= eta )* c(1:mg) )
        if(R_g>0){
          TestStatistic[[g]]$gate2.wg.rej <- ( TestStatistic[[g]]$fdr.j.g <= order.fdr.j.g[R_g] )*1
          TestStatistic[[g]]$eta.g <- sum( order.fdr.j.g[1:R_g] )/R_g
          Rgs[g] <- R_g
        }else{
          TestStatistic[[g]]$gate2.wg.rej <- ( TestStatistic[[g]]$fdr.j.g )*0
          TestStatistic[[g]]$eta.g <- 0
          Rgs[g] <- 0
        }
        TestStatistic[[g]]$fdr.star.g <- 1-( 1- TestStatistic[[g]]$eta.g)*( 1-TestStatistic[[g]]$fdr.g )
        group.fdr.star.g[g] <- TestStatistic[[g]]$fdr.star.g
      }

    ## Step 2. between-group decision
    ORD.g <- order( group.fdr.star.g, decreasing=FALSE )
    number.rej.group <- max ( ( cumsum( Rgs[ORD.g]*group.fdr.star.g[ORD.g] ) <= alpha*cumsum( Rgs[ORD.g] ) ) * c(1:G)  )
    if( number.rej.group>0 )
      {
        ind.rej.group <- ( group.fdr.star.g <= group.fdr.star.g[ ORD.g[ number.rej.group ] ] )*1
        for(g in 1:G){
          if( sum(TestStatistic[[g]]$gate2.wg.rej)!=0){
            TestStatistic[[g]]$gate2.bg.rej <- ind.rej.group[g]
            TestStatistic[[g]]$gate2.wg.rej <- TestStatistic[[g]]$gate2.wg.rej * ind.rej.group[g]
          }else{
            TestStatistic[[g]]$gate2.bg.rej <- 0
            TestStatistic[[g]]$gate2.wg.rej <- TestStatistic[[g]]$gate2.wg.rej * 0
          }
        }
      }else{
        for(g in 1:G){
          TestStatistic[[g]]$gate2.bg.rej <- 0
          TestStatistic[[g]]$gate2.wg.rej <- 0*TestStatistic[[g]]$gate2.wg.rej
        }
      }
    TestStatistic
  }


GATE.3 <- function(TestStatistic, alpha=0.05, eta=alpha)
  {
    
    G <- length(TestStatistic)

    ## A group is rejected if the group local fdr is less than or euqal to eta
    ## Within the group, the fdr level is adjusted according to (alpha- lfdr_g)/(1-lfdr_g)
    for( g in 1:G )
      {
        TestStatistic[[g]]$gate3.bg.rej <- ( TestStatistic[[g]]$fdr.g < eta )
        if( TestStatistic[[g]]$gate3.bg.rej ==0 ){
          TestStatistic[[g]]$gate3.wg.rej <- 0 * TestStatistic[[g]]$X
        }else{
          lfdr.temp <- TestStatistic[[g]]$fdr.j.g
          thresh.temp <- (alpha - TestStatistic[[g]]$fdr.g )/( 1 - TestStatistic[[g]]$fdr.g )
          lfdr.temp.sort <- sort( lfdr.temp, decreasing=FALSE )
          Rg <- max( ( cumsum( lfdr.temp.sort ) <= thresh.temp * c(1: TestStatistic[[g]]$mg)  ) * c(1:TestStatistic[[g]]$mg)  )
          if( Rg == 0){
            TestStatistic[[g]]$gate3.wg.rej <- 0 * TestStatistic[[g]]$X
            TestStatistic[[g]]$gate3.bg.rej <- 0
          }else{
            TestStatistic[[g]]$gate3.wg.rej <- ( TestStatistic[[g]]$fdr.j.g <= ( lfdr.temp.sort[Rg] ) )
          }
        }
      }
    TestStatistic
  }


GATE.4 <- function(TestStatistic, alpha=0.05, eta=alpha, increments=100)
  {
    
      G <- length(TestStatistic)
      group.fdr <- array(0, G)

      ## Collect all between group local fdr
      for( g in 1:G)
          group.fdr[g] <- TestStatistic[[g]]$fdr.g
      ## Select the rejected groups
      group.fdr.sort <- sort(group.fdr,decreasing=FALSE)
      Sg <- max( ( cumsum( group.fdr.sort ) <= eta * c(1: G)  ) * c(1 :G )  )
      if( Sg==0 ){
          for(g in 1:G)
          {
              TestStatistic[[g]]$gate4.bg.rej <- 0
              TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
              TestStatistic[[g]]$gate4.bg.rej.for.bb <- 0
          }
      }else{
          ## Within group rejection
          eta.Sel <- mean( group.fdr.sort[ 1:Sg ] )
          ## thresh.temp <- (alpha-eta.Sel)/(1-eta.Sel )
          max.rej <- array(0, increments)
          PFDR_Sel <- array(0, increments)
          
          for( i in 1:increments)
          {
              Group_Sel <- 0
              thresh.temp <- alpha*i/increments
              cum_Lfdr_g <- 0
              cum_Lfdrgj <- 0
              
              ## Go to the selected group and make final decision
              for(g in 1:G )
              {
                  if( TestStatistic[[g]]$fdr.g <= group.fdr.sort[Sg] ){
                      lfdr.temp <- TestStatistic[[g]]$fdr.j.g
                      lfdr.temp.sort <- sort( lfdr.temp, decreasing=FALSE )
                      Rg <- max( ( cumsum( lfdr.temp.sort ) <= thresh.temp * c(1: TestStatistic[[g]]$mg)  ) * c(1:TestStatistic[[g]]$mg)  )
                      if( Rg > 0){
                          max.rej[i] <- max.rej[i] + Rg
                          Group_Sel <- Group_Sel + 1
                          PFDR_Sel[i] <- PFDR_Sel[i] + sum( lfdr.temp.sort[ 1:Rg ] )/Rg *( 1-TestStatistic[[g]]$fdr.g) + TestStatistic[[g]]$fdr.g
                      }
                  }
              }
              if( Group_Sel>0)
                  PFDR_Sel[i] <- PFDR_Sel[i]/Group_Sel
          }
          ## End here is the code for running within group rejection.
          
          ## Start here, we choose the optimal threshold and rerun the code above.
          max.ind <- max( (PFDR_Sel<= alpha)*c(1:increments) )
          if(max.ind > 0){
              thresh.temp <- alpha*max.ind/increments
              ## Go to the selected group and make final decision
              for(g in 1:G )
              {
                  if( TestStatistic[[g]]$fdr.g > group.fdr.sort[Sg]  ){
                      TestStatistic[[g]]$gate4.bg.rej <- 0
                      TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
                      TestStatistic[[g]]$gate4.bg.rej.for.bb <- 0
                  }else{
                      lfdr.temp <- TestStatistic[[g]]$fdr.j.g
                      lfdr.temp.sort <- sort( lfdr.temp, decreasing=FALSE )
                      Rg <- max( ( cumsum( lfdr.temp.sort ) <= thresh.temp * c(1: TestStatistic[[g]]$mg)  ) * c(1:TestStatistic[[g]]$mg)  )
                      if( Rg == 0){
                          TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
                          TestStatistic[[g]]$gate4.bg.rej <- 0
                      }else{
                          TestStatistic[[g]]$gate4.wg.rej <- ( TestStatistic[[g]]$fdr.j.g <= ( lfdr.temp.sort[Rg] ) )
                          TestStatistic[[g]]$gate4.bg.rej <- 1
                      }
                      TestStatistic[[g]]$gate4.bg.rej.for.bb <- 1
                  }
              }
          }else{
              for(g in 1:G )
              {
                  TestStatistic[[g]]$gate4.bg.rej <- 0
                  TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
                  TestStatistic[[g]]$gate4.bg.rej.for.bb <- 0
              }
          }
          ## End working on the optimal choice of threshold
      }
      TestStatistic
  }


GATE.4.old <- function(TestStatistic, alpha=0.05, eta=alpha)
  {
    
    G <- length(TestStatistic)
    group.fdr <- array(0, G)

    ## Collect all between group local fdr
    for( g in 1:G)
      group.fdr[g] <- TestStatistic[[g]]$fdr.g
    ## Select the rejected groups
    group.fdr.sort <- sort(group.fdr,decreasing=FALSE)
    Sg <- max( ( cumsum( group.fdr.sort ) <= eta * c(1: G)  ) * c(1 :G )  )
    if( Sg==0 ){
      for(g in 1:G)
        {
          TestStatistic[[g]]$gate4.bg.rej <- 0
          TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
          TestStatistic[[g]]$gate4.bg.rej.for.bb <- 0
        }
    }else{
      eta.Sel <- mean( group.fdr.sort[ 1:Sg ] )
      thresh.temp <- (alpha-eta.Sel)/(1-eta.Sel )
      ## Go to the selected group and make final decision
      for(g in 1:G )
        {
          if( TestStatistic[[g]]$fdr.g > group.fdr.sort[Sg]  ){
            TestStatistic[[g]]$gate4.bg.rej <- 0
            TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
            TestStatistic[[g]]$gate4.bg.rej.for.bb <- 0
          }else{
            lfdr.temp <- TestStatistic[[g]]$fdr.j.g
            lfdr.temp.sort <- sort( lfdr.temp, decreasing=FALSE )
            Rg <- max( ( cumsum( lfdr.temp.sort ) <= thresh.temp * c(1: TestStatistic[[g]]$mg)  ) * c(1:TestStatistic[[g]]$mg)  )
            if( Rg == 0){
              TestStatistic[[g]]$gate4.wg.rej <- 0 * TestStatistic[[g]]$X
              TestStatistic[[g]]$gate4.bg.rej <- 0
            }else{
              TestStatistic[[g]]$gate4.wg.rej <- ( TestStatistic[[g]]$fdr.j.g <= ( lfdr.temp.sort[Rg] ) )
              TestStatistic[[g]]$gate4.bg.rej <- 1
            }
            TestStatistic[[g]]$gate4.bg.rej.for.bb <- 1
          }
        }                
    }
    TestStatistic
  }


GBH <- function(TestStatistic, pi1, pi2.1, L, muL, sigmaL, cL, alpha=0.05)
{
    pi0 <- 1-pi1
    pi2.0 <- 1-pi2.1
    G <- length( TestStatistic)
    mgs <- array(0,G)    
    all.weight.p.value <- c()
    
    for( g in 1:G )
    {
        p.value <- 2*pnorm( -abs( TestStatistic[[g]]$X) )        
        p.value <- pi2.0/pi2.1 * p.value
        TestStatistic[[g]]$p.value <- p.value ## Store the weighted p.value
        all.weight.p.value <- c(all.weight.p.value, p.value)
        mgs[g] <- TestStatistic[[g]]$mg
    }
    alpha.omega <- alpha/( 1- pi2.0)
    N <- sum(mgs)
    all.weight.p.value <- sort( all.weight.p.value )
    Rej <- max( ( all.weight.p.value <= ( alpha.omega * c(1:N)/N ) ) * c(1:N)  )
    if( Rej==0){
        for(g in 1:G)
        {
            TestStatistic[[g]]$gbh.bg.rej <- 0
            TestStatistic[[g]]$gbh.wg.rej <- ( TestStatistic[[g]]$p.value)*0
            TestStatistic[[g]]$p.value <- NULL
        }
    }else{
      thresh <- all.weight.p.value[ Rej ]
      for(g in 1:G)
      {
          TestStatistic[[g]]$gbh.wg.rej <- ( TestStatistic[[g]]$p.value <= thresh )
          TestStatistic[[g]]$gbh.bg.rej <- ( sum(TestStatistic[[g]]$gbh.wg.rej) > 0 )
          TestStatistic[[g]]$p.value <- NULL
      }
    }
    TestStatistic
}


SC <- function(TestStatistic, pi1, pi2.1, L, muL, sigmaL, cL, alpha=0.05)
  {
    pi0 <- 1-pi1
    pi2.0 <- 1-pi2.1
    G <- length( TestStatistic )
    all.lfdr <- c()
    
    for( g in 1:G )
      {
        TestStatistic[[g]]$lfdr.sc <- pi2.0*TestStatistic[[g]]$f0x/ TestStatistic[[g]]$fx      
        all.lfdr <- c(all.lfdr, TestStatistic[[g]]$lfdr.sc )
      }
    
    all.lfdr <- sort( all.lfdr, decreasing=FALSE)
    total.Rej <- max ( ( cumsum( all.lfdr ) < c(1:length(all.lfdr))*alpha ) * c(1:length(all.lfdr) ) )
    
    if( total.Rej ==0 ){
      for(g in 1:G)
        {
          TestStatistic[[g]]$sc.bg.rej <- 0
          TestStatistic[[g]]$sc.wg.rej <- ( TestStatistic[[g]]$lfdr.sc)*0
          TestStatistic[[g]]$lfdr.sc <- NULL
        }
    }else{
      thresh <- all.lfdr[ total.Rej ]
      for(g in 1:G)
        {
          TestStatistic[[g]]$sc.wg.rej <- ( TestStatistic[[g]]$lfdr.sc <= thresh )
          TestStatistic[[g]]$sc.bg.rej <- ( sum(TestStatistic[[g]]$sc.wg.rej) > 0 )
          TestStatistic[[g]]$lfdr.sc <- NULL
        }
    }  
    TestStatistic
  }


Pooled <- function(TestStatistic, pi1, pi2.1, L, muL, sigmaL, cL, alpha=0.05)
  {
    pi0 <- 1-pi1
    pi2.0 <- 1-pi2.1
    G <- length( TestStatistic )
    all.pooled.lfdr <- c()
    
    for( g in 1:G )
      {
        TestStatistic[[g]]$lfdr.pooled <- (1-pi1*pi2.1/(1-pi2.0^TestStatistic[[g]]$mg ) ) * TestStatistic[[g]]$f0x/ ( (1-pi1*pi2.1/(1-pi2.0^TestStatistic[[g]]$mg ) )* TestStatistic[[g]]$f0x + pi1*pi2.1/(1-pi2.0^TestStatistic[[g]]$mg ) *TestStatistic[[g]]$f1x )
        all.pooled.lfdr <- c(all.pooled.lfdr, TestStatistic[[g]]$lfdr.pooled )
      }
    
    all.pooled.lfdr <- sort( all.pooled.lfdr, decreasing=FALSE)
    total.Rej <- max ( ( cumsum( all.pooled.lfdr ) < c(1:length(all.pooled.lfdr))*alpha ) * c(1:length(all.pooled.lfdr) ) )
    
    if( total.Rej ==0 ){
      for(g in 1:G)
        {
          TestStatistic[[g]]$pooled.bg.rej <- 0
          TestStatistic[[g]]$pooled.wg.rej <- ( TestStatistic[[g]]$lfdr.pooled)*0
          TestStatistic[[g]]$lfdr.pooled <- NULL
        }
    }else{
      thresh <- all.pooled.lfdr[ total.Rej ]
      for(g in 1:G)
        {
          TestStatistic[[g]]$pooled.wg.rej <- ( TestStatistic[[g]]$lfdr.pooled <= thresh )
          TestStatistic[[g]]$pooled.bg.rej <- ( sum(TestStatistic[[g]]$pooled.wg.rej) > 0 )
          TestStatistic[[g]]$lfdr.pooled <- NULL
        }
    }  
    TestStatistic
  }



BB.Simes <- function(TestStatistic, pi1=0, pi2.1=0, alpha=0.05)
{
    ## BB method where the groups are selected using Simes' method
    if( (pi1!=0) & (pi2.1!=0 ) )
    {
        ## The default option pi1=0, pi2.1=0 indicate that no information is provided for these two proportions.
        
        G <- length(TestStatistic)
        
        ## Calculate the group-wise rejection using Sime's combination
        ## #####
        ## Calculate Sime's combination pvalue first
        Rej.G.ind <- array(0, G)
        Simes.Group <- array(0, G)
        for( g in 1:G ){
            TestStatistic[[g]]$p.value <- 2*pnorm( -abs( TestStatistic[[g]]$X) )
            sort.p.value <- sort( TestStatistic[[g]]$p.value, decreasing=FALSE )
            Simes.Group[g] <- TestStatistic[[g]]$mg * (1-pi1*pi2.1) * min( sort.p.value/c( 1:TestStatistic[[g]]$mg ) )
        }
        ## Decide the rejection based on the group
        sort.Simes.Group <- sort( Simes.Group, decreasing=FALSE )
        Rej.Group <- max( ( sort.Simes.Group <= ( alpha* c(1:G)/G ) )*c(1:G) )
        if( Rej.Group==0 )
        {
            ## If no rejected group, set all decision as false
            for(g in 1:G)
            {
                TestStatistic[[g]]$bbsimes.bg.rej <- 0
                TestStatistic[[g]]$bbsimes.wg.rej <- TestStatistic[[g]]$X *0
                TestStatistic[[g]]$p.value <- NULL
            }
        }else{
            ## If there exist rejected groups
            Rej.G.ind <-  ( Simes.Group <= sort.Simes.Group[ Rej.Group ] ) 
            
            for( g in 1:G )
            {
                if( Rej.G.ind[g]==TRUE){
                    ## Goto the rejected groups
                    temp <- TestStatistic[[g]]$p.value *( 1- pi1 * pi2.1 )
                    mg.tmp <- TestStatistic[[g]]$mg
                    R.tmp <- max( (sort( temp, decreasing=FALSE) <= ( alpha*Rej.Group/G *c(1:mg.tmp)/mg.tmp ) )* c(1:mg.tmp) )
                    if( R.tmp > 0 ){
                        TestStatistic[[g]]$bbsimes.bg.rej = 1
                        TestStatistic[[g]]$bbsimes.wg.rej = (temp <= R.tmp * alpha * Rej.Group/G/mg.tmp )
                    }else{
                        TestStatistic[[g]]$bbsimes.bg.rej = 0
                        TestStatistic[[g]]$bbsimes.wg.rej = temp *0
                    }
                    TestStatistic[[g]]$p.value=NULL
                }else{
                    ## Go to the non-rejected groups
                    TestStatistic[[g]]$bbsimes.bg.rej=0
                    TestStatistic[[g]]$bbsimes.wg.rej= TestStatistic[[g]]$X*0
                    TestStatistic[[g]]$p.value=NULL
                }                
            }
        }
    }
    TestStatistic
}




BB.Oracle <- function(TestStatistic, pi1=0, pi2.1=0, alpha=0.05)
{
    if( (pi1!=0) &(pi2.1!=0) )
    {
        ## In this function, the rejected groups are determined by the one-way gate 2. (Or one-way gate4 in the first version of the paper).
        G <- length(TestStatistic)
        Rej.G.ind <- array(0, G)
        for( g in 1:G ){
            Rej.G.ind[ g ] <- TestStatistic[[g]]$gate4.bg.rej.for.bb
            TestStatistic[[g]]$bb.bg.rej = 0
            TestStatistic[[g]]$bb.wg.rej = array(0, TestStatistic[[g]]$mg )
        }
        
        Rej.G.ind <- which( Rej.G.ind!=0 )
        total.Rej.G <- length(  Rej.G.ind )
        
        if( total.Rej.G > 0 ){
            for( g in Rej.G.ind )
            {
                temp <- 2*pnorm( -abs( TestStatistic[[g]]$X) )*( 1-pi1*pi2.1 )
                mg.tmp <- TestStatistic[[g]]$mg
                R.tmp <- max( (sort( temp, decreasing=FALSE) <= ( alpha*total.Rej.G/G *c(1:mg.tmp)/mg.tmp ) )* c(1:mg.tmp) )
                if( R.tmp > 0 ){
                    TestStatistic[[g]]$bb.bg.rej = 1
                    TestStatistic[[g]]$bb.wg.rej = (temp <= R.tmp * alpha * total.Rej.G/G/mg.tmp )
                }
            }
        }
    }
    TestStatistic
}


BB <- function(TestStatistic, alpha=0.05)
{
    ## In this function, the rejected groups are determined by the one-way gate 2. (Or one-way gate4 in the first version of the paper).
    G <- length(TestStatistic)
    Rej.G.ind <- array(0, G)
    for( g in 1:G ){
        Rej.G.ind[ g ] <- TestStatistic[[g]]$gate4.bg.rej.for.bb
        TestStatistic[[g]]$bb.bg.rej = 0
        TestStatistic[[g]]$bb.wg.rej = array(0, TestStatistic[[g]]$mg )
    }
    
    Rej.G.ind <- which( Rej.G.ind!=0 )
    total.Rej.G <- length(  Rej.G.ind )
    
    if( total.Rej.G > 0 ){
        for( g in Rej.G.ind )
        {
            temp <- 2*pnorm( -abs( TestStatistic[[g]]$X) )
            mg.tmp <- TestStatistic[[g]]$mg
            R.tmp <- max( (sort( temp, decreasing=FALSE) <= ( alpha*total.Rej.G/G *c(1:mg.tmp)/mg.tmp ) )* c(1:mg.tmp) )
            if( R.tmp > 0 ){
                TestStatistic[[g]]$bb.bg.rej = 1
                TestStatistic[[g]]$bb.wg.rej = (temp <= R.tmp * alpha * total.Rej.G/G/mg.tmp )
            }
        }
    }
    TestStatistic
}



GATE.wrapper <- function( TestStatistic, alpha=0.05, eta=alpha, pi1.ini=0.7, pi2.1.ini=0.4, L=2, muL.ini=c(-1,1), sigmaL.ini=c(1,1), cL.ini=c(0.5,0.5), DELTA=0.001, sigma.KNOWN=FALSE )
  {
    em.esti <- GATE.em( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN ) ## Calculate the estimator based on EM algorithm
    TestStatistic <- GATE.localfdr(TestStatistic, em.esti$pi1, em.esti$pi2.1, L, em.esti$muL, em.esti$sigmaL, em.esti$cL) ## Calculate the local fdr scores for all the hypothesis
    GroupTest <- GATE.TLTA(TestStatistic, alpha, eta=alpha)
    GroupTest$parameter <- em.esti
    GroupTest
  }

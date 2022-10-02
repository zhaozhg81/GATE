Assign_Category <- function( X, K, muK, sigmaSqK, mix.prob )
{
    multi.norm.prob <- array(0, K)
    num <- mix.prob * dnorm(X, muK, sqrt(sigmaSqK ) )
    den <- sum( num )
    if(den>0){
        multi.norm.prob <- num/den
        category <- max( rmultinom(1,1,multi.norm.prob)*c(1:K) )
        category
    }else{
        category <- sample(1, c(1:K) )
    }
    category
}

d.Mix.Norm <- function(X, K, muK, sigmaSqK, mix.prob)
{
    sum( mix.prob * dnorm(X, muK, sqrt(sigmaSqK) ) )
}

GATE.Gibbs.Dirichelet <- function( TestStatistic, burnin, thin, size, K=2, d.dirichlet=c(1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE ){

    G <- length( TestStatistic )
    mgs <- array(0, G)
    pi1.post <- array(0, size)
    pi2.post <- array(0, size)
    mu.post <- array(0, c(size, K) )
    sigmaSq.post <- array(0, c(size, K) )
    mix.prob.post <- array(0, c(size, K) )

    ## Initialization
    pi1.tmp <- rbeta(1, alpha1, beta1)
    pi2.tmp <- rbeta(1, alpha2, beta2)
    sigmaSq.tmp <- rgamma(K, 1, 1)
    mu.tmp <- rnorm(K, 1 ,1 )
    mix.prob.tmp <- rdirichlet(1, d.dirichlet )   
    ## mu.tmp <- array(2, K)
    ## sigmaSq.tmp <- array(1, K)
    ## mix.prob.tmp <- c(1)

    
    for(g in 1:G){
        mgs[g] <- TestStatistic[[g]]$mg
    }

    
    for( num.rep in 1:( burnin + thin*size ) )
    ## for(num.rep in 2:100)
    {
        alpha1_prime <- alpha1; beta1_prime <- beta1; alpha2_prime <- alpha2; beta2_prime <- beta2;
            
        ## Step one, update theta.g, update theta.jg
        ## ########################################
        mu.post.mean <- array(0, K)  ## This is the variable for calculating the posterior distribution of mu, \sum_{gj} \theta_{gj}X_{gj}
        total.no.zero <- array(0, K) ## This is the quantity:  \sum_{gj}\theta_{gj}
        sigma.post.scale <- array(0, K)  ## This is the quantity: \sum_{gj}(x_{gj}-mu)^2
        sigma.post.shape <- array(0, K)
        
        ## ## Update theta.g, and theta.jg
        for(g in 1:G)
        {
            ## Get the mixture component
            TestStatistic[[g]]$mix.comp <- array(0, mgs[g] )
            TestStatistic[[g]]$mix.comp <- unlist( lapply( TestStatistic[[g]]$X, Assign_Category, K, mu.tmp, sigmaSq.tmp, mix.prob.tmp ) )

            ## ## Update theta.g
            TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X ) * 10
            ## TestStatistic[[g]]$f1x <- dnorm( TestStatistic[[g]]$X, mu.tmp[ TestStatistic[[g]]$mix.comp], sqrt( sigmaSq.tmp[ TestStatistic[[g]]$mix.comp] ) ) * 10
            TestStatistic[[g]]$f1x <- unlist( lapply( TestStatistic[[g]]$X, d.Mix.Norm, K, mu.tmp, sigmaSq.tmp, mix.prob.tmp ) ) * 10            
            TestStatistic[[g]]$f0prod <- prod( TestStatistic[[g]]$f0x )

            
            post.prob <- (1-pi1.tmp) * TestStatistic[[g]]$f0prod/( (1-pi1.tmp)*TestStatistic[[g]]$f0prod + pi1.tmp* ( prod( (1-pi2.tmp)*TestStatistic[[g]]$f0x + pi2.tmp*TestStatistic[[g]]$f1x) - (1-pi2.tmp)^(mgs[g]) * TestStatistic[[g]]$f0prod )/( 1- (1-pi2.tmp)^(mgs[g]) ) )
            TestStatistic[[g]]$theta.g.tmp <- ( runif(1) > post.prob )          
                  
            ## ## Update theta.jg
            TestStatistic[[g]]$theta.gj.tmp <- array(0, mgs[g] )
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
            
            for(kkk in 1:K){
                mu.post.mean[kkk] <- mu.post.mean[kkk] + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * TestStatistic[[g]]$X * (TestStatistic[[g]]$mix.comp==kkk) ) ## This is the variable for calculating the posterior distribution of mu, \sum_{gj} \theta_{gj}X_{gj}
                total.no.zero[kkk] <- total.no.zero[kkk] + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * (TestStatistic[[g]]$mix.comp==kkk) ) ## This is the quantity:  \sum_{gj}\theta_{gj}
                sigma.post.scale[kkk] <- sigma.post.scale[kkk] + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * ( TestStatistic[[g]]$X-mu.tmp[kkk])^2 * (TestStatistic[[g]]$mix.comp==kkk) ) ## This is the quantity: \sum_{gj}(x_{gj}-mu)^2
                sigma.post.shape[kkk] <- sigma.post.shape[kkk] + TestStatistic[[g]]$theta.g.tmp * sum( TestStatistic[[g]]$theta.gj.tmp * (TestStatistic[[g]]$mix.comp==kkk) )
            }
  
            
            alpha1_prime <- alpha1_prime + TestStatistic[[g]]$theta.g.tmp
            beta1_prime <- beta1_prime + 1 - TestStatistic[[g]]$theta.g.tmp
            alpha2_prime <- alpha2_prime + TestStatistic[[g]]$theta.g.tmp* sum( TestStatistic[[g]]$theta.gj.tmp)
            beta2_prime <- beta2_prime + TestStatistic[[g]]$theta.g.tmp * (mgs[g]-sum(TestStatistic[[g]]$theta.gj.tmp) )
          
        }
        ## Update pi1.tmp and pi2.tmp
        pi1.tmp <- rbeta(1, alpha1_prime, beta1_prime)
        pi2.tmp <- rbeta(1, alpha2_prime, beta2_prime)

    
        ## Update mu.tmp and sigmaSq.tmp
        for(kkk in 1:K)
        {
            if( total.no.zero[kkk] > 0)
            {
                shrinkage.factor <- (1/prec.mu)/( 1/prec.mu + sigmaSq.tmp[kkk]/total.no.zero[kkk])
                mu.tmp[kkk] <- rnorm(1, shrinkage.factor * mu.post.mean[kkk]/total.no.zero[kkk], shrinkage.factor * sigmaSq.tmp[kkk]/total.no.zero[kkk] )
                sigmaSq.tmp[kkk] <- 1/( rgamma(1, shape= shape.sigma + sigma.post.shape[kkk]/2, scale=1/( 1/scale.sigma + sigma.post.scale[kkk]/2 ) ) )
            }
            if( total.no.zero[kkk] == 0)
            {
                mu.tmp[kkk] <- rnorm(1, 0, 1)
                sigmaSq.tmp[kkk] <- 1
            }
        }

        ## Update eta_K for the dirichlet distribution
        mix.prob.tmp <- rdirichlet( 1, d.dirichlet + total.no.zero )


        ## Test
        ## mu.tmp <- muL
        sigmaSq.tmp <- array(1, K)
        ## mix.prob.tmp <- cL
        
        if( ((num.rep-burnin)%%thin ==0 ) && (num.rep>burnin) )
        {
            kk = (num.rep-burnin)/thin
            pi1.post[kk] <- pi1.tmp
            pi2.post[kk] <- pi2.tmp
            mu.post[kk,] <- mu.tmp
            sigmaSq.post[kk,] <- sigmaSq.tmp
            mix.prob.post[kk, ] <- mix.prob.tmp
        }
        if( verbose==TRUE)
        if( num.rep%%1000==0)
            print( paste("Total number of rep: ", burnin+thin*size, "; Current rep: ", num.rep, sep="" ) )

      }
   
    TestStatistic <- GATE.localfdr.Dir(TestStatistic, K, median(pi1.post), median(pi2.post), apply(mu.post, 2, median), apply(sigmaSq.post,2,median), apply( mix.prob.post,2,median )  )
    list( TestStatistic=TestStatistic, mu.post=mu.post, pi1.post=pi1.post, pi2.post=pi2.post, sigmaSq.post=sigmaSq.post, mix.prob.post=mix.prob.post)
    
}


GATE.localfdr.Dir <- function(TestStatistic, K, pi1, pi2.1, mu, sigmaSq, mix.prob)
{
  pi0 <- 1-pi1
  pi2.0 <- 1-pi2.1
  G <- length( TestStatistic )
  for( g in 1:G )
    {
        mg <- TestStatistic[[g]]$mg
        TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X )
        TestStatistic[[g]]$f1x <- unlist( lapply( TestStatistic[[g]]$X, d.Mix.Norm, K, mu, sigmaSq, mix.prob ) ) 
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

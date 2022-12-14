## library(GroupTest)
library(MCMCpack)

source("R/GATE.R")
source("R/Gibbs_Dirichelet.R")



pi1.ini <- 0.5
pi2.1.ini <- 0.6

if(L==1)
{
  muL <- c(2)
  sigmaL <- c(1)
  cL <- c(1)

  muL.ini <- c(0)
  sigmaL.ini <- c(1)
  cL.ini <- c(1)
}

if(L==2)
  {
    muL <- c(-2, 2.5)
    sigmaL <- c( 1.2, 0.8 )
    cL <- c(0.4, 0.6)
    
    muL.ini <- c(1,-1)
    sigmaL.ini <- c(1,1 )
    cL.ini <- c(0.5,0.5)
  }

if(L==3)
  {
    muL <- c(-2, 1, 3.5)
    sigmaL <- c(1.2, 1, 0.8)
    cL <- c(0.4, 0.4, 0.2)
    
    muL.ini <- c(-2,0,2)
    sigmaL.ini <- c(1,1,1)
    cL <- c(0.3,0.3,0.4)
}


DELTA <- 0.001
sigma.KNOWN=FALSE


alpha <- 0.05

data <- array( list(), G )

burnin <- 1000
size <- 500
thin <- 10
numSim <- 100

## burnin <- 10
## size <- 5
## thin <- 2
## numSim <- 100


filename <- paste("result/Gate1_G_",G, "_L_",L, "_mg_",mg, "_pi_BG_",round(pi.BG,digits=3),"_pi_WG_", round(pi.WG, digits=3), ".Rdata",sep="")

FDP <- array(0, c(7, numSim ) )
Rej <- array(0, c(7, numSim ) )
True.Rej <- array(0, c(7, numSim ) )

FDP.G <- array(0, c(7, numSim ) )
Rej.G <- array(0, c(7, numSim ) )
True.Rej.G <- array(0, c(7,numSim ) )
pfdr <- array(0, c(numSim ) )


for(numsim in 1:numSim )
{
    set.seed(numsim)
    
    theta.G <- ( runif( G, 0, 1) < pi.BG )
    X <- array(0, c(G, mg) )
    
    for( g in 1:G )
    {
        data[[g]]$true.theta.G <- theta.G[ g ]
        data[[g]]$mg <- mg
        data[[g]]$true.theta.wg <- ( runif( data[[g]]$mg ) < pi.WG ) * data[[g]]$true.theta.G
        if( theta.G[g] ==0 ){
            data[[g]]$true.theta.wg <- runif( data[[g]]$mg )*0
        }else{
            temp <- ( runif( data[[g]]$mg ) < pi.WG )
            while( sum(temp)==0 )
            {
              temp <- ( runif( data[[g]]$mg ) < pi.WG )
            }
            data[[g]]$true.theta.wg <- temp
        }           
        data[[g]]$X <- (data[[g]]$true.theta.wg == 0 )*rnorm( data[[g]]$mg ) + ( data[[g]]$true.theta.wg==1 )* gen.mix.normal(muL, sigmaL, cL, L, data[[g]]$mg )
        ## data[[g]]$X <- (data[[g]]$true.theta.wg == 0 )*rnorm( data[[g]]$mg ) + ( data[[g]]$true.theta.wg==1 )*rnorm( data[[g]]$mg, muL[1], 1 )
        X[g,] <- data[[g]]$X
    }
    
    data.oracle <- GATE.localfdr(data, pi.BG, pi.WG, L, muL, sigmaL, cL)
    ## data <- GATE.emGibbs( data, burnin, thin, size )

    res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=L, d.dirichlet=rep(1,L), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
    
##    res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=2, d.dirichlet=c(1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
    
    ## res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=3, d.dirichlet=c(1,1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )

    pi1.esti <- median( res.dir$pi1.post )
    pi2.esti <- median( res.dir$pi2.post )
    mu.esti <- apply( res.dir$mu,2, median)
    sigmaSq.esti <- apply( res.dir$sigmaSq.post, 2, median )
    mix.prob.esti <- apply( res.dir$mix.prob.post,2, median )

    

    data <- res.dir$TestStatistic

    ## res <- GATE.Gibbs( data, burnin, thin, size, alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=FALSE )
    ## data <- res$TestStatistic

       
    data.oracle <- GATE.1( data.oracle, alpha=0.05 )
    data <- GATE.1(data, alpha=0.05 )
    data <- SC(data, pi1.esti, pi2.esti, L, mu.esti, sigmaSq.esti, mix.prob.esti, alpha=0.05)
    data <- Pooled(data, pi1.esti, pi2.esti, L, mu.esti, sigmaSq.esti, mix.prob.esti, alpha=0.05)
    data <- GBH(data, pi.BG, pi.WG, L, muL, sigmaL, cL, alpha=0.05)
    data <- GBH.datadriven(data, method="LSL", alpha=0.05)
    data <- GBH.datadriven(data, method="TST", alpha=0.05)
    ## EM Algorithm
    pi2.1.ini <- 0.5
    esti <- GATE.em( data, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN)
    data.em <- GATE.localfdr(data, esti$pi1, esti$pi2.1, L, esti$muL, esti$sigmaL, esti$cL)
    data.em <- GATE.1( data.em, alpha=0.05)
    
    
    for(g in 1:G){
      ## GATE 1
        True.Rej.G[1, numsim ] <- True.Rej.G[1, numsim ] + (theta.G[g]* data[[g]]$gate1.bg.rej  )
        Rej.G[1, numsim] <- Rej.G[1, numsim ] +   data[[g]]$gate1.bg.rej
        True.Rej[1, numsim ] <- True.Rej[1, numsim ] + sum( data[[g]]$true.theta.wg * data[[g]]$gate1.wg.rej )
        Rej[1,  numsim ] <- Rej[1, numsim ] + sum( data[[g]]$gate1.wg.rej )
        
        ## GBH, LSL
        True.Rej.G[2, numsim ] <- True.Rej.G[2, numsim ] + (theta.G[g]* data[[g]]$gbh.lsl.bg.rej  )
        Rej.G[2, numsim] <- Rej.G[2, numsim ] +   data[[g]]$gbh.lsl.bg.rej
        True.Rej[2, numsim ] <- True.Rej[2, numsim ] + sum( data[[g]]$true.theta.wg * data[[g]]$gbh.lsl.wg.rej )
        Rej[2,  numsim ] <- Rej[2, numsim ] + sum( data[[g]]$gbh.lsl.wg.rej )
        
        ## SC
        True.Rej.G[3, numsim ] <- True.Rej.G[3, numsim ] + (theta.G[g]* data[[g]]$sc.bg.rej  )
        Rej.G[3, numsim] <- Rej.G[3, numsim ] +   data[[g]]$sc.bg.rej
        True.Rej[3, numsim ] <- True.Rej[3, numsim ] + sum( data[[g]]$true.theta.wg * data[[g]]$sc.wg.rej )
        Rej[3,  numsim ] <- Rej[3, numsim ] + sum( data[[g]]$sc.wg.rej )
        
        ## Pooled
        True.Rej.G[4, numsim ] <- True.Rej.G[4, numsim ] + (theta.G[g]* data[[g]]$pooled.bg.rej  )
        Rej.G[4, numsim] <- Rej.G[4, numsim ] +   data[[g]]$pooled.bg.rej
        True.Rej[4, numsim ] <- True.Rej[4, numsim ] + sum( data[[g]]$true.theta.wg * data[[g]]$pooled.wg.rej )
        Rej[4,  numsim ] <- Rej[4, numsim ] + sum( data[[g]]$pooled.wg.rej )

        ## GBH,TST
        True.Rej.G[5, numsim ] <- True.Rej.G[5, numsim ] + (theta.G[g]* data[[g]]$gbh.tst.bg.rej  )
        Rej.G[5, numsim] <- Rej.G[5, numsim ] +   data[[g]]$gbh.tst.bg.rej
        True.Rej[5, numsim ] <- True.Rej[5, numsim ] + sum( data[[g]]$true.theta.wg * data[[g]]$gbh.tst.wg.rej )
        Rej[5,  numsim ] <- Rej[5, numsim ] + sum( data[[g]]$gbh.tst.wg.rej )

        ## GATE Oracle
        True.Rej.G[6, numsim ] <- True.Rej.G[6, numsim ] + (theta.G[g]* data.oracle[[g]]$gate1.bg.rej  )
        Rej.G[6, numsim] <- Rej.G[6, numsim ] +   data.oracle[[g]]$gate1.bg.rej
        True.Rej[6, numsim ] <- True.Rej[6, numsim ] + sum( data.oracle[[g]]$true.theta.wg * data.oracle[[g]]$gate1.wg.rej )
        Rej[6,  numsim ] <- Rej[6, numsim ] + sum( data.oracle[[g]]$gate1.wg.rej )

        ## EM algorithm
        True.Rej.G[7, numsim ] <- True.Rej.G[7, numsim ] + (theta.G[g]* data.em[[g]]$gate1.bg.rej  )
        Rej.G[7, numsim] <- Rej.G[7, numsim ] +   data.em[[g]]$gate1.bg.rej
        True.Rej[7, numsim ] <- True.Rej[7, numsim ] + sum( data.em[[g]]$true.theta.wg * data.em[[g]]$gate1.wg.rej )
        Rej[7,  numsim ] <- Rej[7, numsim ] + sum( data.em[[g]]$gate1.wg.rej )

        
    }
    print( paste("lambda_ind=", lambda.ind, "; numsim=", numsim, sep="" ) )
}

save.image(file=filename)

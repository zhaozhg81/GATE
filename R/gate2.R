## library(GroupTest)
library(MCMCpack)

source("R/GATE.R")
source("R/Gibbs_Dirichelet.R")

numSim <- 100
DELTA <- 0.001

alpha <- 0.05
data <- array( list(), G )

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


eta.all <- c(1:16)/16*alpha
alpha_star_all <- array(0, c(4,numSim, length(eta.all) ) )


FDP <- array(0, c(5, numSim, length(eta.all) ) )
Rej <- array(0, c(5, numSim, length(eta.all) ) )
True.Rej <- array(0, c(5, numSim, length(eta.all) ) )

FDP.G <- array(0, c(5, numSim, length(eta.all) ) )
Rej.G <- array(0, c(5, numSim, length(eta.all) ) )
True.Rej.G <- array(0, c(5,numSim, length(eta.all) ) )
pfdr <- array(0, c(numSim, length(eta.all) ) )

SelFDR <- array( 0, c(5,numSim, length(eta.all) ) )
Ave.SelFDR <- array(0, c(5,numSim, length(eta.all) ) )

burnin <- 1000
size <- 500
thin <- 10
numSim <- 100

filename <- paste("result/Gate2_G_",G, "_mg_",mg[1],"_L_",L, "_pi_BG_",pi.BG,"_pi_WG_",pi.WG, ".Rdata",sep="")

##for(numsim in (5*(ind-1)+1):(5*ind) )
for(numsim in 1:numSim )
{
    theta.G <- ( runif( G, 0, 1) < pi.BG )

    
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
        data[[g]]$X <- (data[[g]]$true.theta.wg == 0 )*rnorm( data[[g]]$mg ) + ( data[[g]]$true.theta.wg==1 )*gen.mix.normal(muL, sigmaL, cL, L, data[[g]]$mg )
    }
    
    data.oracle <- GATE.localfdr(data, pi.BG, pi.WG, L, muL, sigmaL, cL)
    
##     res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=1, d.dirichlet=c(1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
    res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=L, d.dirichlet=rep(1, L), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
    ## res.dir <- GATE.Gibbs.Dirichelet( data, burnin, thin, size, K=3, d.dirichlet=c(1,1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
    data <- res.dir$TestStatistic

    
    ## res.embayes <- GATE.emGibbs( data, burnin, thin, size, alpha1=1, beta1=1, alpha2=1, beta2=1, verbose=FALSE )
    data.embayes <- res.dir$TestStatistic


    ## EM Algorithm
    esti <- GATE.em( data, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA)
    data.em <- GATE.localfdr(data, esti$pi1, esti$pi2.1, L, esti$muL, esti$sigmaL, esti$cL)



    for(eta.ind in c( 1:length(eta.all) ) )
    {
        eta <- eta.all[ eta.ind ]
        
        onewaygate4 <- GATE.4(data, alpha=0.05, eta=eta )
        data <- onewaygate4$TestStatistic
        alpha_star_all[1, numsim, eta.ind] = onewaygate4$alpha_star
          
        onewaygate4.oracle <- GATE.4(data.oracle,alpha=0.05,eta=eta)
        data.oracle <- onewaygate4.oracle$TestStatistic
        alpha_star_all[2,numsim, eta.ind] = onewaygate4.oracle$alpha_star
        
        data <- BB.Simes(data, pi.BG, pi.WG, alpha=0.05)
        data <- BB.Oracle(data, pi.BG, pi.WG, alpha=0.05)
        gate4.embayes <-GATE.4(data.embayes,alpha=0.05,eta=eta)
        data.embayes <- gate4.embayes$TestStatistic
        alpha_star_all[3,numsim,eta.ind] <- gate4.embayes$alpha_star
        
        gate4.em <- GATE.4(data.em,alpha=0.05,eta=eta)
        data.em <- gate4.em$TestStatistic
        alpha_star_all[4,numsim,eta.ind] <- gate4.em$alpha_star
        
        Rg.gate4 <- 0
        Rg.oracle.gate4 <- 0
        Rg.bbsimes <- 0
        Rg.bb <- 0
        Rg.embayes <- 0
        Rg.em <- 0

        
        for(g in 1:G){
            True.Rej.G[1, numsim, eta.ind ] <- True.Rej.G[1, numsim, eta.ind ] + (theta.G[g] * data[[g]]$gate4.bg.rej  )
            Rej.G[1, numsim, eta.ind] <- Rej.G[1, numsim, eta.ind ] +   data[[g]]$gate4.bg.rej
            True.Rej[1, numsim, eta.ind ] <- True.Rej[1, numsim, eta.ind ] + sum( data[[g]]$true.theta.wg * data[[g]]$gate4.wg.rej )
            Rej[1,  numsim, eta.ind ] <- Rej[1, numsim, eta.ind ] + sum( data[[g]]$gate4.wg.rej )
            SelFDR[1, numsim, eta.ind ] <- SelFDR[ 1,numsim, eta.ind ] + sum( (1-data[[g]]$true.theta.wg) * data[[g]]$gate4.wg.rej )/max( sum( data[[g]]$gate4.wg.rej), 1)
            Rg.gate4 <- Rg.gate4 + data[[g]]$gate4.bg.rej

            
            True.Rej.G[2, numsim, eta.ind ] <- True.Rej.G[2, numsim, eta.ind ] + (theta.G[g] * data[[g]]$bb.bg.rej  )
            Rej.G[2, numsim, eta.ind] <- Rej.G[2, numsim, eta.ind ] +   data[[g]]$bb.bg.rej
            True.Rej[2, numsim, eta.ind ] <- True.Rej[2, numsim, eta.ind ] + sum( data[[g]]$true.theta.wg * data[[g]]$bb.wg.rej )
            Rej[2,  numsim, eta.ind ] <- Rej[2, numsim, eta.ind ] + sum( data[[g]]$bb.wg.rej )
            SelFDR[2, numsim, eta.ind ] <- SelFDR[ 2,numsim, eta.ind ] + sum( (1-data[[g]]$true.theta.wg) * data[[g]]$bb.wg.rej )/max( sum( data[[g]]$bb.wg.rej), 1)
            Rg.bb <- Rg.bb + data[[g]]$bb.bg.rej
            
            True.Rej.G[3, numsim, eta.ind ] <- True.Rej.G[3, numsim, eta.ind ] + (theta.G[g] * data[[g]]$bbsimes.bg.rej  )
            Rej.G[3, numsim, eta.ind] <- Rej.G[3, numsim, eta.ind ] +   data[[g]]$bbsimes.bg.rej
            True.Rej[3, numsim, eta.ind ] <- True.Rej[3, numsim, eta.ind ] + sum( data[[g]]$true.theta.wg * data[[g]]$bbsimes.wg.rej )
            Rej[3,  numsim, eta.ind ] <- Rej[3, numsim, eta.ind ] + sum( data[[g]]$bbsimes.wg.rej )
            SelFDR[3, numsim, eta.ind ] <- SelFDR[ 3,numsim, eta.ind ] + sum( (1-data[[g]]$true.theta.wg) * data[[g]]$bbsimes.wg.rej )/max( sum( data[[g]]$bbsimes.wg.rej), 1)
            Rg.bbsimes <- Rg.bbsimes + data[[g]]$bbsimes.bg.rej

            True.Rej.G[4, numsim, eta.ind ] <- True.Rej.G[4, numsim, eta.ind ] + (theta.G[g] * data.oracle[[g]]$gate4.bg.rej  )
            Rej.G[4, numsim, eta.ind] <- Rej.G[4, numsim, eta.ind ] +   data.oracle[[g]]$gate4.bg.rej
            True.Rej[4, numsim, eta.ind ] <- True.Rej[4, numsim, eta.ind ] + sum( data.oracle[[g]]$true.theta.wg * data.oracle[[g]]$gate4.wg.rej )
            Rej[4,  numsim, eta.ind ] <- Rej[4, numsim, eta.ind ] + sum( data.oracle[[g]]$gate4.wg.rej )
            SelFDR[4, numsim, eta.ind ] <- SelFDR[ 4,numsim, eta.ind ] + sum( (1-data.oracle[[g]]$true.theta.wg) * data.oracle[[g]]$gate4.wg.rej )/max( sum( data.oracle[[g]]$gate4.wg.rej), 1)
            Rg.oracle.gate4 <- Rg.oracle.gate4 + data.oracle[[g]]$gate4.bg.rej



            True.Rej.G[5, numsim, eta.ind ] <- True.Rej.G[5, numsim, eta.ind ] + (theta.G[g] * data.em[[g]]$gate4.bg.rej  )
            Rej.G[5, numsim, eta.ind] <- Rej.G[5, numsim, eta.ind ] +   data.em[[g]]$gate4.bg.rej
            True.Rej[5, numsim, eta.ind ] <- True.Rej[5, numsim, eta.ind ] + sum( data.em[[g]]$true.theta.wg * data.em[[g]]$gate4.wg.rej )
            Rej[5,  numsim, eta.ind ] <- Rej[5, numsim, eta.ind ] + sum( data.em[[g]]$gate4.wg.rej )
            SelFDR[5, numsim, eta.ind ] <- SelFDR[ 5,numsim, eta.ind ] + sum( (1-data.em[[g]]$true.theta.wg) * data.em[[g]]$gate4.wg.rej )/max( sum( data.em[[g]]$gate4.wg.rej), 1)
            Rg.em <- Rg.em + data.em[[g]]$gate4.bg.rej            

        }
        
        SelFDR[1,numsim,eta.ind] <- SelFDR[1,numsim, eta.ind]/Rg.gate4       
        SelFDR[2,numsim,eta.ind] <- SelFDR[2,numsim, eta.ind]/Rg.bb
        SelFDR[3,numsim,eta.ind] <- SelFDR[3,numsim, eta.ind]/Rg.bbsimes
        SelFDR[4,numsim,eta.ind] <- SelFDR[4,numsim, eta.ind]/Rg.oracle.gate4
        SelFDR[5,numsim,eta.ind] <- SelFDR[5,numsim, eta.ind]/Rg.em
        
        print( paste("eta=",eta, "; numsim=", numsim, sep="" ) )
    }
}

save.image(file=filename)

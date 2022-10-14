library(MCMCpack)
library(coda)
library(lattice)
library(ggplot2)
library(gridExtra)

data(AYP)

source("../R/GATE.R")
source("../R/Gibbs_Dirichelet.R")

fig_dir =  "../figure/realdata"



gen.mix.normal <- function(muL, sigmaL, cL, L, mg)
{
    x <- array(0, mg)
    for( i in 1:mg ){
        xL <- rnorm( L, muL, sigmaL )
        x[i] <- xL[ sum( runif(1) > cumsum( cL ) )  + 1 ]
    }
    x
  }

dmixnormal <- function(X, L, muL, sigmaL, mix.prob)
  {
    n <- length(X)
    den <- array(0, n)
    for(i in 1:n)
      {
        den[i] <- sum( mix.prob * dnorm( X[i], muL, sqrt(sigmaL) ) )
      }
    den
  }


bh.cut <- function(pvalue, alpha)
{
    p <- length(pvalue)
    R <- max(  ( sort(pvalue, decreasing=FALSE) <= (alpha*c(1:p)/p) )*c(1:p) )
    if(R==0)
    {
        thresh <- -1
    }else{
        thresh <- R*alpha/p
    }
    thresh
}


alpha <- 0.05
G <- length(AYP)

x.all <- c()
for(g in 1:G)
  x.all <- append(x.all, AYP[[g]]$X)

if( !file.exists("AYP_Bayes_Analysis.RData") )
{
    
    burnin <- 10000
    size <- 500
    thin <- 20
    
    if( !file.exists("AYP_Chains.Rdata") ){
        ## Generate three MCMC chains 
        res1.dir <- GATE.Gibbs.Dirichelet( AYP, burnin, thin, size, K=2, d.dirichlet=c(1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
        res2.dir <- GATE.Gibbs.Dirichelet( AYP, burnin, thin, size, K=2, d.dirichlet=c(1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
        res3.dir <- GATE.Gibbs.Dirichelet( AYP, burnin, thin, size, K=2, d.dirichlet=c(1,1), alpha1=1, beta1=1, alpha2=1, beta2=1, prec.mu=0.0001, shape.sigma=0.0001, scale.sigma=1/0.0001, verbose=TRUE )
        save.image("AYP_Chains.Rdata")
    }else{
        load("AYP_Chains.Rdata")
      }

    
    res1.dir$mu.post <- res1.dir$mu.post[, c(2,1) ]
    res1.dir$mix.prob.post <- res1.dir$mix.prob.post[, c(2,1) ]
    
    ## Check the mixing
    x1 <- cbind( res1.dir$pi1.post, res1.dir$pi2.post, res1.dir$mu.post, res1.dir$mix.prob )
    colnames(x1) <- c("pi1","pi2","mu1","mu2","c1","c2")
    mc.x1 <- mcmc(x1)
    
    x2 <- cbind( res2.dir$pi1.post, res2.dir$pi2.post, res2.dir$mu.post, res2.dir$mix.prob )
    colnames(x2) <- c("pi1","pi2","mu1","mu2","c1","c2")
    mc.x2 <- mcmc(x2)
    
    x3 <- cbind( res3.dir$pi1.post, res3.dir$pi2.post, res3.dir$mu.post, res3.dir$mix.prob )
    colnames(x3) <- c("pi1","pi2","mu1","mu2","c1","c2")
    mc.x3 <- mcmc(x3)
    
    mclist.AYP <- as.mcmc.list(list(mc.x1,mc.x2,mc.x3) )
    densityplot( mclist.AYP )
    
    
    
    ## Run GATE 1
    
    
    ## Get the estimated parameter
    pi1 <- median( c( res1.dir$pi1.post, res2.dir$pi1.post, res3.dir$pi1.post ) )
    pi2 <- median( c( res1.dir$pi2.post, res2.dir$pi2.post, res3.dir$pi2.post ) )
    L <- 2
    muL <- apply( rbind( res1.dir$mu.post, res2.dir$mu.post, res3.dir$mu.post), 2, median )
    sigmaL <- c(1, 1)
    mix.prob <- apply( rbind( res1.dir$mix.prob, res2.dir$mix.prob, res3.dir$mu.post), 2, median )
        
    AYP <- GATE.localfdr(AYP, pi1, pi2, L, muL, sigmaL, mix.prob)
    AYP <- GATE.1(AYP, alpha=0.05 )
    AYP <- GATE.4(AYP, alpha=0.05 )
    
    Rej.G <- 0
    Rej <- 0

    Rej.4 <- 0
    Rej.G.4 <- 0
    
    for( g in 1:G)
    {
        Rej.G <- Rej.G + AYP[[g]]$gate1.bg.rej
        Rej <- Rej + sum( AYP[[g]]$gate1.wg.rej )

        Rej.G.4 <- Rej.G.4+ AYP[[g]]$gate4.bg.rej
        Rej.4 <- Rej.4 + sum( AYP[[g]]$gate4.wg.rej )
    }

    save.image("AYP_Bayes_Analysis.RData")
}else{
  load("AYP_Bayes_Analysis.RData")
}


## Check the model assumption
pi1 <- median( res1.dir$pi1.post)
pi2 <- median( res1.dir$pi2.post)
L <- 1
muL <- apply( res1.dir$mu.post,2,median)
sigmaL <- c(1,1)
mix.prob <- apply( res1.dir$mix.prob.post,2,median)

xx <- c(1:2000)/100-10
yy <- (1-pi1)*dnorm(xx) + pi1*dmixnormal(xx, L, muL, sigmaL, mix.prob)

postscript("fittedCurve.eps", horizontal=FALSE)

LWD <- 2
CEX <- 2
plot( density( x.all), ylim=c(0,0.3), col='black', lwd=LWD, main="Density", xlab="x", ylab="pdf", cex.main=CEX)
points(xx,yy, 'l', col='red', lwd=LWD)
legend(x=1, y=0.3, c("Kernel Density", "Estimated"), lty=c(1,1), col=c('black','red'), lwd=c(LWD,LWD) )

dev.off()




total = 0
for(g in 1:G)
  {
    if( !is.na(pmatch("New Haven", AYP[[g]]$School.District)  ) )
      print( paste("New Haven School District:", g, sep="") )

    if( !is.na(pmatch("Berkeley Unified", AYP[[g]]$School.District)  ) )
      print( paste("Berkeley Unified School District:", g, sep="") )
    
  }

L = 2
AYP <- SC(AYP, pi1, pi2, L, muL, sigmaL, mix.prob, alpha=0.05)
## AYP <- SC(AYP, pi1, pi2, L, muL, c(1,1), mix.prob , alpha=0.05)
## data <- GBH.datadriven(AYP, method="LSL", alpha=0.05)
AYP <- GBH.datadriven(AYP, method="TST", alpha=0.05)

## AYP <- GBH(AYP, pi1, pi2, L, muL, c(1,1), mix.prob, alpha=0.05)
AYP <- Pooled(AYP, pi1, pi2, L, muL, c(1,1), mix.prob, alpha=0.05)


L=2
pi2.1.ini <- 0.5
pi1.ini <- 0.5

  muL <- c(-2, 2.5)
  sigmaL <- c( 1.2, 0.8 )
  cL <- c(0.4, 0.6)
  
  muL.ini <- c(1,-1)
  sigmaL.ini <- c(1,1 )
  cL.ini <- c(0.5,0.5)

  DELTA=0.001
sigma.KNOWN=TRUE 

esti <- GATE.em( data, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN)
AYP.em <- GATE.localfdr(AYP, esti$pi1, esti$pi2.1, L, esti$muL, esti$sigmaL, esti$cL)
AYP.em <- GATE.1( AYP.em, alpha=0.05)


rej.sc <- 0
rej.gbh <- 0
rej.pooled <- 0
rej.gate.em <- 0
z.all <- c()
decision.gate.all <- c()
decision.sc.all <- c()
decision.gbh.all <- c()
decision.pooled.all <- c()
decision.gate.em.all <- c()


for(g in 1:G)
  {
    rej.sc <- rej.sc + AYP[[g]]$sc.bg.rej * sum( AYP[[g]]$sc.wg.rej )
    rej.gate.em <- rej.gate.em + AYP.em[[g]]$gate1.bg.rej * sum( AYP.em[[g]]$gate1.wg.rej)
    rej.gbh <- rej.gbh + AYP[[g]]$gbh.tst.bg.rej * sum( AYP[[g]]$gbh.tst.wg.rej )
    rej.pooled <- rej.pooled + AYP[[g]]$pooled.bg.rej * sum( AYP[[g]]$pooled.wg.rej )

    z.all <- append(z.all, AYP[[g]]$X )
    decision.gate.all <- append(decision.gate.all, AYP[[g]]$gate1.bg.rej*AYP[[g]]$gate1.wg.rej )
    decision.gate.em.all <- append(decision.gate.em.all, AYP.em[[g]]$gate1.bg.rej*AYP[[g]]$gate1.wg.rej )
    decision.sc.all <- append(decision.sc.all, AYP[[g]]$sc.bg.rej * AYP[[g]]$sc.wg.rej )
    decision.gbh.all <- append( decision.gbh.all, AYP[[g]]$gbh.tst.bg.rej * AYP[[g]]$gbh.tst.wg.rej)
    decision.pooled.all <- append(decision.pooled.all, AYP[[g]]$pooled.bg.rej * AYP[[g]]$pooled.wg.rej )
  }

## source("processdata.R")


p <- 4118
combined <- read.csv("combined_processed.csv",sep=",",header=T)
combined$gate <- decision.gate.all
combined$sc <- decision.sc.all
combined$gbh <- decision.gbh.all
combined$pooled <- decision.pooled.all
combined$gate.em <- decision.gate.em.all


exclude.ind <- union( which( is.na(combined$z) ), which(is.na(combined$z.1) ) )
exclude.ind <- union( exclude.ind, which( abs(combined$z-z.all)>0.001) )
exclude.ind <- union( exclude.ind,  which( abs(combined$z) < 2 )  )
exclude.ind <- union( exclude.ind, which( abs(combined$z) > 5 ) )
exclude.ind <- union( exclude.ind, which( abs(combined$z.1) > 5 ) )

keep.ind <- setdiff( c(1:p), exclude.ind )

DEC <- c('rej', 'fail-to-rejct')

combined.keep <- combined[keep.ind,]
combined.keep$gate <- DEC[ 2- combined.keep$gate ]
combined.keep$gate.em <- DEC[ 2- combined.keep$gate.em ]
combined.keep$sc <- DEC[ 2- combined.keep$sc ]
combined.keep$gbh <- DEC[ 2- combined.keep$gbh ]
combined.keep$pooled <- DEC[ 2- combined.keep$pooled ]


postscript(paste(fig_dir,"AYP_gate.eps",sep=""),  horizontal=FALSE)
gplot.gate <- ggplot( combined.keep, aes(z, z.1, colour=gate,shape=gate ) ) + geom_point()
gplot.gate <- gplot.gate + geom_hline(yintercept=2, linetype="dashed", color='black')
gplot.gate <- gplot.gate + geom_hline(yintercept=-2, linetype="dashed", color='black')
gplot.gate <- gplot.gate  + labs(x="AYP 2013") + labs(y="AYP 2015")
gplot.gate
dev.off()

postscript(paste(fig_dir, "AYP_gate_em.eps",sep=""),  horizontal=FALSE)
gplot.gate <- ggplot( combined.keep, aes(z, z.1, colour=gate.em,shape=gate.em ) ) + geom_point()
gplot.gate <- gplot.gate + geom_hline(yintercept=2, linetype="dashed", color='black')
gplot.gate <- gplot.gate + geom_hline(yintercept=-2, linetype="dashed", color='black')
gplot.gate <- gplot.gate  + labs(x="AYP 2013") + labs(y="AYP 2015")
gplot.gate
dev.off()


postscript( paste(fig_dir,"AYP_sc.eps",sep=""),  horizontal=FALSE)
gplot.sc <- ggplot( combined.keep, aes(z, z.1, colour=sc,shape=sc ) ) + geom_point()
gplot.sc <- gplot.sc + geom_hline(yintercept=2, linetype="dashed", color='black')
gplot.sc <- gplot.sc + geom_hline(yintercept=-2, linetype="dashed", color='black')
gplot.sc <- gplot.sc  + labs(x="AYP 2013") + labs(y="AYP 2015")
gplot.sc
dev.off()


postscript(paste(fig_dir,"AYP_pooled.eps",sep=""),  horizontal=FALSE)
gplot.pooled <- ggplot( combined.keep, aes(z, z.1, colour=pooled,shape=pooled ) ) + geom_point()
gplot.pooled <- gplot.pooled + geom_hline(yintercept=2, linetype="dashed", color='black')
gplot.pooled <- gplot.pooled + geom_hline(yintercept=-2, linetype="dashed", color='black')
gplot.pooled <- gplot.pooled  + labs(x="AYP 2013") + labs(y="AYP 2015")
gplot.pooled
dev.off()

## combined.keep$naive <- combined.keep$pooled
## postscript("AYP_naive.eps",  horizontal=FALSE)
## gplot.naive <- ggplot( combined.keep, aes(z, z.1, colour=naive,shape=naive ) ) + geom_point()
## gplot.naive <- gplot.naive + geom_hline(yintercept=2, linetype="dashed", color='black')
## gplot.naive <- gplot.naive + geom_hline(yintercept=-2, linetype="dashed", color='black')
## gplot.naive <- gplot.naive  + labs(x="AYP 2013") + labs(y="AYP 2015")
## gplot.naive
## dev.off()



postscript(paste(fig_dir,"AYP_gbh.eps",sep=""), horizontal=FALSE)
gplot.gbh <- ggplot( combined.keep, aes(z, z.1, colour=gbh,shape=gbh ) ) + geom_point()
gplot.gbh <- gplot.gbh + geom_hline(yintercept=2, linetype="dashed", color='black')
gplot.gbh <- gplot.gbh + geom_hline(yintercept=-2, linetype="dashed", color='black')
gplot.gbh <- gplot.gbh + labs(x="AYP 2013") + labs(y="AYP 2015")
gplot.gbh
dev.off()


rej.ind.gate = which( combined.keep$gate=='rej')

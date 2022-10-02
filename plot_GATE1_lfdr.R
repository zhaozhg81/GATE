figdir <- "./"

library(latex2exp)

numSim <- 100

G <- 100
mg <- 50
L <- 1

pi.WGs = c(0.1,0.3,0.5,0.7)
Gs = c(100, 1000)
mgs = c(50, 5)
Ls =c(1,2)

for( ii in 1:4)
  for(jj in 1:2)
    for(kk in 1:2)
    {
      L = Ls[kk]
      G = Gs[jj]
      mg = mgs[jj]

      lfdr.error <- array(0, c(9) )
      total.count <- array(0, c(9) )
      
      
      for(lambda.IND in 1:9 )
        {
          total.count[ lambda.IND ] <- 0
          
          pi.BGs <- c(1:9)/10
          pi.BG <- pi.BGs[lambda.IND]
          
          if( mg<=10 )
            {
              
              lambdas <- c( 1:9)^2/49
              lambda <- lambdas[ lambda.IND ]    
              pi.WG <- pi.WGs[ii]
              
              tmp <- (1-(1-pi.WG)^mg[1])/( (1-pi.WG)^mg[1]) * lambda
              pi.BG <- tmp/(1+tmp)
              
            }else{     
              pi.BGs <- c(1:9)/10      
              pi.WG <- pi.WGs[ii]
              pi.BG <- pi.BGs[lambda.IND]
            }
          
          filename <- paste("result/Gate1_G_",G,"_L_",L, "_mg_",mg,"_pi_BG_",round(pi.BG,digits=3),"_pi_WG_", round(pi.WG, digits=3), ".Rdata",sep="")
          load(filename)
          
          
          temp <- GATE.localfdr.Bayes( res.dir$TestStatistic, pi.BG, pi.WG, muL, sigmaL)
          
          for(g in 1:G)
            {          
              lfdr.error[ lambda.IND ] <- lfdr.error[lambda.IND] +   sum( (temp[[g]]$fdr.j.g<0.1)* ( res.dir$TestStatistic[[g]]$fdr.j.g - temp[[g]]$fdr.j.g)^2 )
              total.count[ lambda.IND ] <-  total.count[lambda.IND] + sum( temp[[g]]$fdr.j.g<0.1 )
            }
        }
    
      
      mean.lfdr.error = apply( lfdr.error/total.count, 1, mean )
      
      ## Method lists
      ## 1. Gate 1
      ## 2. GBH, LSL
      ## 3. SC
      ## 4. Naive
      ## 5. GBH, TST
      ## 6. GATE Oracle
      ## 7. GATE EM Algorithm
      
      IND=c(1:9)
      if( (G== 1000) & (mg==5) & (L==2) & (pi.WG==0.7) )
        IND = c(1:4, 6:9)
      if( (G== 100) & (mg==50) & (L==2) & (pi.WG==0.3) )
        IND = c(1:8)
      if( (G== 100) & (mg==50) & (L==2) & (pi.WG==0.5) )
        IND = c(1:7,9)
      
      
      mu <- 2
      
      CEX <- 1.5
      LWD <- 3
      
      if( mg[1] <= 10 ){
        
        
        lfdr.figure <- paste(figdir,"figure/Bayes_lfdr_gate1_G_",G,"_mg_",mg[1], "_L_",L, "_pi_2_",pi.WG,".eps",sep="")
        
        postscript( file=lfdr.figure, horizontal=FALSE)
        
        
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas[IND], mean.lfdr.error[IND], type="b", pch=2, col='green', main=TeX(paste("Squard Error of $\\widehat{Lfdr}$: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'),ylab="Squared Error", ylim=c(0, 1.1* max( mean.lfdr.error) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        
        
        dev.off()
        
      }else{
        
        
        lfdr.figure <- paste(figdir,"figure/Bayes_lfdr_gate1_G_",G,"_mg_",mg[1], "_L_",L, "_pi_2_",pi.WG,".eps",sep="")
        postscript( file=lfdr.figure, horizontal=FALSE)
        
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs[IND], mean.lfdr.error[IND], type="b", pch=2, col='green', main=TeX(paste("Squard Error of $\\widehat{Lfdr}$: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="Squared Error", ylim=c(0, 1.1* max( mean.lfdr.error) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        
        
        dev.off()
      }
    }

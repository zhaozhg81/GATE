figdir <- "~/Dropbox/Apps/Overleaf/GATE_EJS/"

figdir <- "./"

library(latex2exp)

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
      
      numSim <- 100
      
      FDR <- array(0, c(7, 9) )
      True.R <- FDR
      FDR.G <- FDR
      True.R.G <- FDR
      R <- FDR
      R.G <- FDR   
      
      for(lambda.IND in 1:9 )
        {
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
          
          
          filename <- paste("result/Gate1_G_",G, "_L_",L, "_mg_",mg, "_pi_BG_",round(pi.BG,digits=3),"_pi_WG_", round(pi.WG, digits=3), ".Rdata",sep="")
          
          load(filename)
          
          FDR[, lambda.IND] <- apply( (Rej-True.Rej)/(Rej+(Rej==0)),1,mean)
          True.R[,lambda.IND] <- apply( True.Rej,1, mean)
          FDR.G[, lambda.IND] <- 1-apply(True.Rej.G/( Rej.G+ (Rej.G==0) ),1,mean)
          True.R.G[,lambda.IND] <- apply( True.Rej.G,1,mean)
          R[,lambda.IND] <- apply( Rej,1,mean )
          R.G[, lambda.IND ] <- apply( True.Rej.G,1,mean)  
        }
      
      ## Method lists
      ## 1. Gate 1
      ## 2. GBH, LSL
      ## 3. SC
      ## 4. Naive
      ## 5. GBH, TST
      ## 6. GATE Oracle
      ## 7. GATE EM Algorithm
      
      FDR[ is.na(FDR) ] <- 0
      True.R[ is.na(True.R) ] <- 0
      FDR.G[ is.na(FDR.G) ] <- 0
      True.R.G[ is.na(True.R.G) ] <- 0
      R[ is.na(R) ] <- 0
      R.G[ is.na(R.G) ] <- 0
      mu <- 2
      
      CEX <- 1.5
      LWD <- 3
      
      if( mg[1] <= 10 ){
        
        fdr.figure <- paste(figdir,"figure/Bayes_fdr_gate1_G_",G,"_mg_",mg[1],"_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        g.fdr.figure <- paste(figdir,"figure/Bayes_group_fdr_gate1_G_",G,"_mg_",mg[1], "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        
        true.rej.figure <- paste(figdir,"figure/Bayes_true_rej_gate1_G_",G,"_mg_",mg[1],  "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        true.g.rej.figure <- paste(figdir,"figure/Bayes_true_group_rej_gate1_G_",G,"_mg_",mg[1], "_L_",L, "_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        
        rej.figure <- paste(figdir,"figure/Bayes_rej_gate1_G_",G,"_mg_",mg[1],  "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        g.rej.figure <- paste(figdir,"figure/Bayes_group_rej_gate1_G_",G,"_mg_",mg[1], "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
          
        
        postscript( file=fdr.figure, horizontal=FALSE)

        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, FDR[5, ], type="b", pch=2, col='green', main=TeX(paste("FDR: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'),ylab="FDR", ylim=c(0, 1.1* max( FDR) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, FDR[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, FDR[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        points( lambdas, 0*lambdas+ alpha, type="b", pch=25, col='black', lwd=LWD )
        ##        points( lambdas, FDR[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, FDR[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, FDR[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, FDR[1,], type="b", pch=22, col='red', lwd=LWD)
        
        
        legend( "topright", y=max(FDR),  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        
        dev.off()
        
        postscript( file=g.fdr.figure, horizontal=FALSE)
        
    
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, FDR.G[5, ], type="b", pch=2, col='green', main=TeX(paste("FDRG: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'), ylab="Group FDR", ylim=c(0, 1.1* max( FDR.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, FDR.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, FDR.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        points( lambdas, 0*lambdas+ alpha, type="b", pch=25, col='black', lwd=LWD )
        ##        points( lambdas, FDR.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, FDR.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, FDR.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, FDR.G[1,], type="b", pch=22, col='red', lwd=LWD)
        
        
        legend( "topright", y=max(FDR.G),  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        
        dev.off()
        
        
        postscript( file=rej.figure, horizontal=FALSE)
        
        
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, R[5, ], type="b", pch=2, col='green', main=TeX(paste("Rejection: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'), ylab="Rejection", ylim=c(0, max( R) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, R[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, R[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( lambdas, R[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, R[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, R[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, R[1,], type="b", pch=22, col='red', lwd=LWD)
        
        ## Method lists
        ## 1. Gate 1
        ## 2. GBH, LSL
        ## 3. SC
        ## 4. Naive
        ## 5. GBH, TST
        ## 6. GATE Oracle
        ## 7. GATE EM Algorithm
        
        legend( "topleft", y=max(R)/2, c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        
        
        dev.off()
        
    
        postscript( file=true.rej.figure, horizontal=FALSE)
        
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, True.R[5, ], type="b", pch=2, col='green', main=TeX(paste("True Rejection: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'), ylab="True Rej", ylim=c(0, max( True.R) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, True.R[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, True.R[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( lambdas, True.R[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, True.R[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, True.R[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, True.R[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend( "topleft", y=max(True.R)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        dev.off()
        
        
        postscript( file=g.rej.figure, horizontal=FALSE)
        
    
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, R.G[5, ], type="b", pch=2, col='green', main=TeX(paste("Rejected Groups: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'), ylab="Total Group Rej", ylim=c(0, max( R.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, R.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, R.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( lambdas, R.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, R.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, R.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, R.G[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend( "topleft", y=max(R.G)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        
        
        dev.off()
        
        postscript( file=true.g.rej.figure, horizontal=FALSE)
        
    
        
        par(mar=c(5,5,4.1,1.5))
        plot( lambdas, True.R.G[5, ], type="b", pch=2, col='green', main=TeX(paste("Truly Rejected Groups: $\\pi_2=", pi.WG) ), xlab=TeX('$\\lambda$'), ylab="True Group Rej", ylim=c(0, max( True.R.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( lambdas, True.R.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( lambdas, True.R.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( lambdas, True.R.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( lambdas, True.R.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( lambdas, True.R.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( lambdas, True.R.G[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend( "topleft", y=max(True.R.G)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )

        

        dev.off()
        
    
      }else{

        fdr.figure <- paste(figdir,"figure/Bayes_fdr_gate1_G_",G,"_mg_",mg[1],"_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        g.fdr.figure <- paste(figdir,"figure/Bayes_group_fdr_gate1_G_",G,"_mg_",mg[1], "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        
        true.rej.figure <- paste(figdir,"figure/Bayes_true_rej_gate1_G_",G,"_mg_",mg[1],  "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        true.g.rej.figure <- paste(figdir,"figure/Bayes_true_group_rej_gate1_G_",G,"_mg_",mg[1], "_L_",L, "_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        
        rej.figure <- paste(figdir,"figure/Bayes_rej_gate1_G_",G,"_mg_",mg[1],  "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")
        g.rej.figure <- paste(figdir,"figure/Bayes_group_rej_gate1_G_",G,"_mg_",mg[1], "_L_",L,"_pi_WG_", round(pi.WG, digits=3),".eps",sep="")                      
  
        postscript( file=fdr.figure, horizontal=FALSE)
        

        ## Method lists
        ## 1. Gate 1
        ## 2. GBH, LSL
        ## 3. SC
        ## 4. Naive
        ## 5. GBH, TST
        ## 6. GATE Oracle
        ## 7. GATE EM Algorithm
        
    
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, FDR[5, ], type="b", pch=2, col='green', main=TeX(paste("FDR: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="FDR", ylim=c(0, 1.1* max( FDR) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, FDR[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, FDR[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        points( pi.BGs, 0*pi.BGs+ alpha, type="b", pch=25, col='black', lwd=LWD )
        ## points( pi.BGs, FDR[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, FDR[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, FDR[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, FDR[1,], type="b", pch=22, col='red', lwd=LWD)
        

        legend( "topright", y=max(FDR),  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )

        dev.off()
        
        postscript( file=g.fdr.figure, horizontal=FALSE)
        
        
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, FDR.G[5, ], type="b", pch=2, col='green', main=TeX(paste("Group FDR: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="Group FDR", ylim=c(0, 1.1* max( FDR.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, FDR.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, FDR.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        points( pi.BGs, 0*pi.BGs+ alpha, type="b", pch=25, col='black', lwd=LWD )
        ## points( pi.BGs, FDR.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, FDR.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, FDR.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, FDR.G[1,], type="b", pch=22, col='red', lwd=LWD)
        
        
        legend("topright", y=max(FDR.G),  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )

        dev.off()
        

        postscript( file=rej.figure, horizontal=FALSE)
        
    
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, R[5, ], type="b", pch=2, col='green', main=TeX(paste("Rejection: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'),ylab="Rejection", ylim=c(0, max( R) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, R[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, R[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( pi.BGs, R[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, R[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, R[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, R[1,], type="b", pch=22, col='red', lwd=LWD)
        
        
        legend("topleft", y=max(R)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )

        
        dev.off()

    
        postscript( file=true.rej.figure, horizontal=FALSE)
        
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, True.R[5, ], type="b", pch=2, col='green', main=TeX(paste("True Rejection: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="True Rej", ylim=c(0, max( True.R) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, True.R[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, True.R[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( pi.BGs, True.R[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, True.R[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, True.R[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, True.R[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend("topleft", y=max(True.R)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        dev.off()
        
        
        postscript( file=g.rej.figure, horizontal=FALSE)

        
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, R.G[5, ], type="b", pch=2, col='green', main=TeX(paste("Rejected Groups: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="Total Group Rej", ylim=c(0, max( R.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, R.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, R.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ##  points( pi.BGs, R.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, R.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, R.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, R.G[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend("topleft", y=max(R.G)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )
        
        
        dev.off()
        
        postscript( file=true.g.rej.figure, horizontal=FALSE)   
        
        par(mar=c(5,5,4.1,1.5))
        plot( pi.BGs, True.R.G[5, ], type="b", pch=2, col='green', main=TeX(paste("Truly Rejected Groups: $\\pi_2=", pi.WG) ), xlab=TeX('$\\pi_1$'), ylab="True Group Rej", ylim=c(0, max( True.R.G) ), lwd=LWD, cex.lab=CEX, cex.main=CEX )
        points( pi.BGs, True.R.G[3, ], type="b", pch=17, col='blue', lwd=LWD )
        points( pi.BGs, True.R.G[4, ], type="b", pch=0, col='mediumorchid', lwd=LWD )
        ## points( pi.BGs, True.R.G[5,], type="b", pch=21, col='skyblue', lwd=LWD)
        points( pi.BGs, True.R.G[7,], type="b", pch=8, col='darkviolet', lwd=LWD)
        points( pi.BGs, True.R.G[6,], type="b", pch=10, col='coral2', lwd=LWD)
        points( pi.BGs, True.R.G[1,], type="b", pch=22, col='red', lwd=LWD)    
        legend("topleft", y=max(True.R.G)/2,  c("GBH","SC","Naive","GATE, EM", "GATE, Oracle", "GATE"), pch=c(2, 17,0,8,10,22), col=c("green","blue","mediumorchid","darkviolet","coral2","red"), cex=CEX, lwd=LWD )



        dev.off()
        
      }
    }

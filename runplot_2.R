figdir <- "~/Dropbox/Apps/Overleaf/GATE_EJS/figure/"
## figdir <- "./figure/"

numSim <- 100

FDR.all <- array(0, c(5, 3, 16) )
True.R.all <- FDR.all
FDR.G.all <- FDR.all
True.R.G.all <- FDR.all
R.all <- FDR.all
R.G.all <- FDR.all
SelFDR.all <- FDR.all
FDP.all <- FDR.all



L <- 2

G <- 1000
mg <- array(5,  G)
PIWG <- 0.3
alpha_star = array(0, c(3,4,100,16) )

for( lambda.IND in c(1:3) )  
{
  pi.BGs <- 3*c(1:3)/10
  
  pi.WG <- PIWG
  pi.BG <- pi.BGs[lambda.IND]
  
  filename <- paste("result/Gate2_G_",G, "_mg_",mg[1],"_L_",L, "_pi_BG_",pi.BG,"_pi_WG_",pi.WG, ".Rdata",sep="")
          
  
  load(filename)
  alpha_star[lambda.IND,,,] = alpha_star_all[ ,,  ]
  
  for( iii in 1:16)
    {
      True.R.all[,lambda.IND,iii] <- True.R.all[, lambda.IND, iii] + apply( True.Rej[, ,iii], 1, sum)/100
      R.all[,lambda.IND,iii] <- R.all[, lambda.IND, iii] + apply( Rej[, ,iii], 1, sum)/100
      R.G.all[,lambda.IND,iii] <- R.G.all[, lambda.IND, iii] + apply( Rej.G[, ,iii], 1, sum)/100
      True.R.G.all[,lambda.IND,iii] <- True.R.G.all[, lambda.IND, iii] + apply( True.Rej.G[, ,iii], 1, sum)/100
      SelFDR.all[,lambda.IND,iii] <- SelFDR.all[, lambda.IND, iii] + apply( SelFDR[, ,iii], 1, sum)/100
      FDR.all[,lambda.IND,iii] <- FDR.all[,lambda.IND,iii] + apply( ( Rej[, ,iii] - True.Rej[, ,iii])/(Rej[, ,iii] + (Rej[, ,iii]==0) ), 1, sum )/100
      FDR.G.all[,lambda.IND,iii] <- FDR.G.all[,lambda.IND,iii] + apply( ( Rej.G[, ,iii] - True.Rej.G[, ,iii])/(Rej.G[, ,iii] + (Rej.G[, ,iii]==0) ), 1, sum )/100
    }
}

alpha_star_gate2 = array(0, c(3,16) )
alpha_star_gate2_em = array(0, c(3,16)) 

for( iii in 1:16)
{
  alpha_star_gate2[,iii] = apply(alpha_star[,1,,iii],1,mean)
  alpha_star_gate2_em[,iii] = apply(alpha_star[,4,,iii],1,mean)
  
}


mu.pos.all <- c(15:30)/10
CEX <- 2
LWD <- 3




for( iii in 1:length( pi.BGs ) )
{
    if( mg[1]<=10 )
    {
        
        lambdas <- c( 1:20)^2/100
        
        pi.WG <- PIWG
        
        tmp <- (1-(1-pi.WG)^mg[1])/( (1-pi.WG)^mg[1]) * lambdas
        pi.BGs <- tmp/(1+tmp)

        pi.BG <- pi.BGs[iii]
        oracle=TRUE
        
      }else{
          
        lambdas <- c( 1:20)^2/100
        lambda <- lambdas[ iii ]    
        pi.BGs <- c(1:20)/21
        
        pi.WG <- PIWG
        pi.BG <- pi.BGs[ iii ]
        oracle=TRUE
    }
    

    fdr.figure <- paste(figdir,"Bayes_fdr_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),"_L_",L, ".eps",sep="")
    g.fdr.figure <- paste(figdir,"Bayes_group_fdr_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG, digits=3), "_piWG_", 1000*round(pi.WG,digits=3), "_L_",L,  ".eps",sep="")
    sel.fdr.figure <- paste(figdir,"Bayes_sel_fdr_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG, digits=3), "_piWG_", 1000*round(pi.WG,digits=3), "_L_",L, ".eps",sep="")
    
    true.rej.figure <- paste(figdir,"Bayes_true_rej_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),  "_L_",L, ".eps",sep="")
    true.g.rej.figure <- paste(figdir,"Bayes_true_group_rej_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),  "_L_",L, ".eps",sep="")
    
    rej.figure <- paste(figdir,"Bayes_rej_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),  "_L_",L,  ".eps",sep="")
    g.rej.figure <- paste(figdir,"Bayes_group_rej_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),  "_L_",L,  ".eps",sep="")

    
    alpha.star.figure <- paste(figdir,"alpha_star_gate2_G_",G,"_mg_",mg[1],"_piBG_", 1000*round(pi.BG,digits=3), "_piWG_", 1000*round(pi.WG,digits=3),  "_L_",L,  ".eps",sep="")

    postscript( file=alpha.star.figure, horizontal=FALSE)
    
    plot( eta.all, alpha_star_gate2[iii,], 'l', col='red', xlab=expression(eta), ylim=c(0, 0.1), ylab="alpha_star", cex=CEX, lwd=LWD )
    points( eta.all, R.G.all[1,iii,]/G *alpha , 'l', col='black', cex=CEX, lwd=LWD)
    legend("topright", c("alpha star","alpha*Rg/G"), 
           pch=c(1,1), col=c("red","black"), cex=CEX, lwd=LWD )
    
    dev.off()
    
            
    postscript( file=fdr.figure, horizontal=FALSE)
    
    par(mar=c(4,5,4,0))
    
    plot( eta.all, FDR.all[1, iii, ], type="b", pch=22, col='red', main="FDR", xlab=expression(eta), ylab="FDR", ylim=c(0, max( FDR.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, FDR.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, FDR.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, FDR.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, FDR.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    points( eta.all, 0.05* array(1, length(eta.all)), type="b", col='black', pch=1, lwd=LWD)
##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)","alpha"), 
       ##    pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet","black"), cex=CEX, lwd=LWD )
    
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)","alpha"), 
           pch=c(22,0,8,1), col=c("red","blue","darkviolet","black"), cex=CEX, lwd=LWD )
    
    dev.off()

    postscript( file=sel.fdr.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))
    
    plot( eta.all, SelFDR.all[1, iii, ], type="b", pch=22, col='red', main="Sel FDR", xlab=expression(eta), ylab="Sel FDR", 
          ylim=c(0, max( SelFDR.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
  ##  points( eta.all, SelFDR.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, SelFDR.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
  ##  points( eta.all, SelFDR.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, SelFDR.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    points( eta.all, 0.05* array(1, length(eta.all)), type="b", col='black', pch=1, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","Gate 2 (EM)","alpha"), 
    ##       pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet","black"), cex=CEX, lwd=LWD )
    
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)","alpha"), 
           pch=c(22,0,8,1), col=c("red","blue","darkviolet","black"), cex=CEX, lwd=LWD )
    
    dev.off()

    
    postscript( file=g.fdr.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))

    plot( eta.all, FDR.G.all[1, iii, ], type="b", pch=22, col='red', main="Group FDR", xlab=expression(eta), ylab="Group FDR", 
          ylim=c(0, max( FDR.G.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, FDR.G.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, FDR.G.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, FDR.G.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, FDR.G.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    points( eta.all, 0.05* array(1, length(eta.all)), type="b", col='black', pch=1, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)","alpha"), 
    ##        pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet","black"), cex=CEX, lwd=LWD )
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)","alpha"), 
           pch=c(22,0,8,1), col=c("red","blue","darkviolet","black"), cex=CEX, lwd=LWD )
    

    dev.off()
    
    
    postscript( file=rej.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))

    plot( eta.all, R.all[1, iii, ], type="b", pch=22, col='red', main="Number of Rejection", xlab=expression(eta), ylab="Rej", 
          ylim=c(0, max( R.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, R.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, R.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, R.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, R.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)"), 
    ##       pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet"), cex=CEX, lwd=LWD )
    
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)"), 
           pch=c(22,0,8), col=c("red","blue","darkviolet"), cex=CEX, lwd=LWD )
    
    dev.off()
    
    postscript( file=true.rej.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))

    plot( eta.all, True.R.all[1, iii, ], type="b", pch=22, col='red', main="Number of True Rejection", xlab=expression(eta), ylab="True Rej", 
          ylim=c(0, max( True.R.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, True.R.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, True.R.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, True.R.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, True.R.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
  ##    legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)"), 
    ##       pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet"), cex=CEX, lwd=LWD )

    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)"), 
           pch=c(22,0,8), col=c("red","blue","darkviolet"), cex=CEX, lwd=LWD )
    
    
    dev.off()
    
    postscript( file=g.rej.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))

    plot( eta.all, R.G.all[1, iii, ], type="b", pch=22, col='red', main="Number of Rejected Group", xlab=expression(eta), ylab="Rej Group", 
          ylim=c(0, max( R.G.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, R.G.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, R.G.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, R.G.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, R.G.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)"), 
    ##       pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet"), cex=CEX, lwd=LWD )
    
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)"), 
           pch=c(22,0,8), col=c("red","blue","darkviolet"), cex=CEX, lwd=LWD )
    
    
    dev.off()
    
    postscript( file=true.g.rej.figure, horizontal=FALSE)
    par(mar=c(4,5,4,0))

    plot( eta.all, True.R.G.all[1, iii, ], type="b", pch=22, col='red', main="Number of Truly Rejected Group", xlab=expression(eta), ylab="Rej True Group", 
          ylim=c(0, max( True.R.G.all[,iii,], 1.5*alpha) ), lwd=LWD, cex.lab=CEX, cex.main=CEX, lty=1 )
    ## points( eta.all, True.R.G.all[2, iii, ], type="b", pch=17, col='green', lwd=LWD ) 
    points( eta.all, True.R.G.all[3, iii, ], type="b", pch=0,  col='blue', lwd=LWD ) 
    ## points( eta.all, True.R.G.all[4, iii, ], type="b", col='skyblue', pch=10, lwd=LWD ) 
    points( eta.all, True.R.G.all[5, iii, ], type="b", col='darkviolet', pch=8, lwd=LWD)
    ##    points( eta.all, 0*eta.all+ alpha + 1.645 * sqrt( alpha*(1-alpha)/numSim) , 'l', col='black', lwd=LWD, lty=3 )
    ##    points( eta.all, FDR[4, iii,], 'l', col='red', lty=2, lwd=LWD)
    ##    legend(x=0.02, y=0.6*max(FDR[,iii,]), c("Gate","BB", expression(alpha + 1.645* sqrt(alpha*(1-alpha)/100)), "Gate, Oracle"), lty=c(1,2,3,2), col=c("red","green","black","red"), cex=CEX, lwd=LWD )
    ## legend("bottomright", c("GATE 2","BB Simes","BB","GATE 2 (Oracle)","GATE 2 (EM)"), 
    ##       pch=c(22,17,0,10,8,1), col=c("red","green","blue","skyblue","darkviolet"), cex=CEX, lwd=LWD )
    
    legend("bottomright", c("GATE 2","BB","GATE 2 (EM)"), 
           pch=c(22,0,8), col=c("red","blue","darkviolet"), cex=CEX, lwd=LWD )
    

    dev.off()
    
    
}
    

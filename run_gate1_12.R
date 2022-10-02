ind <- 12

para.setting <- data.frame( m = c( rep(100,9), rep(1000,9)), n=c( rep(50,9),rep(5,9)), lambda.ind=rep(c(1:9),2) )

G <- para.setting$m[ind]
mg <- para.setting$n[ind]
lambda.ind <- para.setting$lambda.ind[ind]

for(L in 1:3)
  {

    if( mg <= 10 )
      {        
        lambdas <- c( 1:9)^2/49
        lambda <- lambdas[ lambda.ind ]    
        
        pi.WG <- 0.7
        pi.BG <- 0.3
        
        tmp <- (1-(1-pi.WG)^mg[1])/( (1-pi.WG)^mg[1]) * lambda
        pi.BG <- tmp/(1+tmp)
        
      }else{        
        pi.BGs <- c(1:9)/10
        
        pi.WG <- 0.7
        pi.BG <- pi.BGs[lambda.ind]
        
      }
    source("R/gate1.R")
  }

ind <- 10


para.setting <- data.frame( m = c( rep(500,9), rep(1000,9)), n=c( rep(20,9),rep(5,9)), pi1= rep(c(1,2,3),6), K= rep( c(rep(1,3),rep(2,3),rep(3,3)), 2 ) ) 


##if( mg[1]< 10 )
##{
##   
##    lambdas <- c( 1:20)^2/100
##    lambda <- lambdas[ ind ]    
##    
##    pi.WG <- 0.3
##    pi.BG <- 0.3
##    
##    tmp <- (1-(1-pi.WG)^mg[1])/( (1-pi.WG)^mg[1]) * lambda
##    pi.BG <- tmp/(1+tmp)
##    
##    oracle=FALSE
##    source("gate4.R")
##
##}else{

G <- para.setting$m[ind]
mg <- para.setting$n[ind]


pi.BGs <- 3*c(1:3)/10

pi.WG <- 0.3
pi.BG <- pi.BGs[ para.setting$pi1[ind] ]
L <- para.setting$K[ind]

oracle=FALSE
source("R/gate2.R")

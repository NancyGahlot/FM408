
HestonCallMonteCarlo <-
  function(lambda, vbar, eta, rho, v0, r, tau, S0, K, nSteps=250, nPaths=10000) {{
    
    n <- nSteps
    N <- nPaths
    
    dt <- tau / n
    
    S <- rep(S0,N)
    v <- rep(v0,N)
    V <- c()
    AV <- c()
    AVdev <- c()
    
    for (i in 1:n)
    {
     W1 <- rnorm(N);
     W2 <- rnorm(N);
      W2 <- rho*W1 + sqrt(1 - rho^2)*W2;
      
      sqvdt <- sqrt(v*dt)
      S <- S*exp((r-v/2)*dt + sqrt(v * dt) * W1)
      
      
      
        sqvdt <- sqrt(v*dt)
        v <- v + lambda*(vbar - v)* dt + eta * sqvdt * W2
        v <- ifelse(v<0, -v, v) #to deal with negative variance
    }
    
    
    ## Evaluate mean call value for each path
     
    V <- exp(-r*tau)*(S>K)*(S - K); # Boundary condition for European call
    AV <- mean(V);
    AVdev <- 2 * sd(V) / sqrt(N) 
    
    return(AV)
    #lower = AV-AVdev
    #upper = AV+AVdev
       }
  }


Call_price = rep(0,1)
for (i in 0:1){
  Call_price[i] = HestonCallMonteCarlo(lambda=2, vbar=0.01, eta=0.1, rho=-0.7, v0=0.01, r=0.0,
                                       tau=1, S0=100, K=100)
}
Call_price




BSvol <- function(pc, S, k, price, d, r, t, start= 0.2)
{
  #pc    put/call indicator call=1, put=-1
  #S     Stock price at 0
  #K     strike
  #price option premium
  #d     dividend yield
  #r     riskless rate
  #t     time to maturity
  #start starting value for vol, optional, by default=0.2
  
  
  voli = start
  pricei = BSprice(pc, S, k, voli, d, r, t)
  vegai = BSvega(pc, S, k, voli, d, r, t)
  while(abs(price - pricei) > 0.000001) 
  {
    voli<-voli + (price - pricei) / vegai
    pricei<-BSprice(pc, S, k, voli, d, r, t)
    vegai<-BSvega(pc, S, k, voli, d, r, t)
  }
  IV = BSvol 
  return(BSvol) }

IV <- c()
IV = BSvol(1,100 ,100 ,Call_price, 0, 0, 1/250 )
IV 


BSprice<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  d2 = d1 - vol * sqrt(t)
  
  BSprice = pc * exp(-d * t) * S * 
    pnorm(pc * d1) - pc * k * exp(-r * t) * pnorm(pc * d2)
  return(BSprice)
}


BSdelta<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  
  if (pc == 1) {BSdelta = exp(-d * t) * pnorm(d1)} else 
  {BSdelta = exp(-d * t) * (pnorm(d1) - 1)}
  return(BSdelta)
}


BSvega<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  
  BSvega = exp(-d * t) * S * sqrt(t) * exp((-d1 ^ 2) / 2) / (sqrt(2 * pi))
  return(BSvega)
}


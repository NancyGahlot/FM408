nPaths <- 1000
nSteps_d <- 250
nSteps_h <- 2500
tau <- 1


#########STOCK PRICE############

Stockprice <- function(K,S0,dt,nSteps_d,nPaths, tau){
   
  n <- nSteps_d
  N <- nPaths
  
  dt <- tau / n
  
  S <- rep(S0,N)
   set.seed(1)
   W1 <- rnorm(N)
   
   sqvdt <- sqrt(v*dt)
   S <- S*exp((r-v/2)*dt + sqrt(v * dt) * W1)
  
   return(S)
}

Stock_price_daily <- c()
Stock_price_daily = Stockprice(100,100,dt,nSteps_d, nPaths,1)
Stock_price_daily

######HESTON#########

HestonCallMonteCarlo <-
  function(lambda, vbar, eta, rho, v0, r, tau, S0, K, nSteps_d= 250, nPaths=1000) {
    
    n <- nSteps_d
    N <- nPaths
    
    dt <- tau / n
    
    S <- rep(S0,N)
    v <- rep(v0,N)
    V <- c()
    AV <- c()
    AVdev <- c()
    
  
    set.seed(1)
      W1 <- rnorm(N);
      set.seed(2)
      W2 <- rnorm(N);
      W2 <- rho*W1 + sqrt(1 - rho^2)*W2;
      
      sqvdt <- sqrt(v*dt)
      S <- S*exp((r-v/2)*dt + sqrt(v * dt) * W1)
      
      
      sqvdt <- sqrt(v*dt)
      v <- v + lambda*(vbar - v)* dt + eta * sqvdt * W2
      v <- ifelse(v<0, -v, v) #to deal with negative variance
    
    
    ##Evaluate mean call value for each path
    
    V <- exp(-r*tau)*(S>K)*(S - K); # Boundary condition for European call
    AV <- mean(V);
    AVdev <- 2 * sd(V) / sqrt(N) 
    
    return(V)
    }
  


Call_price_daily = c()
Call_price_daily = HestonCallMonteCarlo(lambda=2, vbar=0.01, eta=0.1, rho= -0.7, v0=0.01, r=0.0,
                                       tau=1, S0=100, K=100)

Call_price_daily




##### BS code ########

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


BSprice(1,100,100,IV_daily,0,0,1)


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


BSvega(1,100,100,IV_daily,0,0,1)


BSvol<-function(pc, S, k, price, d, r, t, start = 0.03){
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
    while(abs(price - pricei) > 0.00001) 
    {
      voli<-voli + (price - pricei) / vegai
      pricei<-BSprice(pc, S, k, voli, d, r, t)
      vegai<-BSvega(pc, S, k, voli, d, r, t)
    }
  }
  BSvol = voli
  return(BSvol)
}

IV_daily <-c()
IV_daily = BSvol(1,Stock_price_daily,100,Call_price_daily, 0, 0, 1, 0.03)
IV_daily




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

BSdelta(1,100,100,IV_daily,0,0,1)



BSgamma<-function(pc, S, k, vol, d, r, t)
{
  
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  
  BSgamma = exp(-d * t) * exp((-d1 ^ 2) / 2) / (sqrt(2 * pi) * S * vol * sqrt(t))
  
  return(BSgamma)}

BSgamma(1,100,100,IV_daily,0,0,1)



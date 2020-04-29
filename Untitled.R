
## BS code

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




BSvol<-function(pc, S, k, price, d, r, t, start = 0.2)
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
  
  BSvol = voli
  return(BSvol)
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

BStheta<-function(pc, S, k, vol, d, r, t) #there is a q in the formula?
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
  
  BStheta = -exp(-d * t) * exp((-d1 ^ 2) / 2) * S * vol / 
    (sqrt(2 * pi) * 2 * sqrt(t)) + pc * d * S * exp(-d * t) * pnorm(pc * d1) - pc * r * k * exp(-r * t) * pnorm(pc * d2)
  return(BStheta)
}

BSrho<-function(pc, S, k, vol, d, r, t)
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
  
  BSrho = pc * k * t * exp(-r * t) * pnorm(pc * d2)
  return(BSrho)
}


###HESTON###

# Risk-Neutral Probability P1

HestonP1<-function(phi, kappa, theta, lambda, rho, sigma, tau, K, S, r, v)
{
  mu1 = 0.5;
  b1 = (kappa + lambda - rho * sigma);
  d1 = sqrt((complex(real=0,imaginary=(rho*sigma*phi)) - b1)^2 - (complex(real=0,imaginary=(sigma^2 * 2 * mu1 * phi)) - sigma^2 * phi^2));
  g1 = (b1 - complex(real=0,imaginary=(rho * sigma * phi)) + d1) / (b1 - complex(real=0,imaginary=(rho * sigma * phi)) - d1);
  DD1_1 = (b1 - complex(real=0,imaginary=(rho * sigma * phi)) + d1) / (sigma^2);
  DD1_2 = 1 - exp(d1 * tau);
  DD1_3 = 1 - g1 * exp(d1 * tau);
  DD1 = DD1_1 * (DD1_2 / DD1_3);
  CC1_1 = complex(real=0,imaginary=(r * phi * tau));
  CC1_2 = (kappa * theta) / (sigma^2);
  CC1_3 = (b1 - complex(real=0,imaginary=(rho * sigma * phi)) + d1) * tau;
  CC1_4 = 2 * log((1 - g1 * exp(d1 * tau)) / (1 - g1));
  cc1 = CC1_1 + CC1_2 * (CC1_3 - CC1_4);
  f1 = exp(cc1 + DD1 * v + complex(real=0,imaginary=phi * log(S)));
  
  y = Re(f1 * exp(complex(real=0,imaginary=-phi * log(K))) / (complex(real=0,imaginary=phi)));
  return(y)
}




# Risk-Neutral Probability P2

HestonP2<-function(phi, kappa, theta, lambda, rho, sigma, tau, K, S, r, v)
{
  mu1 = -0.5;
  b1 = kappa + lambda;
  d1 = sqrt((rho * sigma * complex(real=0,imaginary=phi) - b1)^2 - (sigma^2 * 2 * mu1 * complex(real=0,imaginary=phi ) - sigma^2 * phi^2));
  g1 = (b1 - rho * sigma * complex(real=0,imaginary=phi ) + d1) / (b1 - rho * sigma * complex(real=0,imaginary=phi ) - d1);
  
  DD1_1 = (b1 - rho * sigma * complex(real=0,imaginary=phi ) + d1) / (sigma^2);
  DD1_2 = 1 - exp(d1 * tau);
  DD1_3 = 1 - g1 * exp(d1 * tau);
  DD1 = DD1_1 * DD1_2 / DD1_3;
  
  CC1_1 = complex(real=0,imaginary=r * phi * tau);
  CC1_2 = kappa * theta / (sigma^2);
  CC1_3 = (b1 - rho * sigma * complex(real=0,imaginary=phi ) + d1) * tau;
  CC1_4 = 2 * log((1 - g1 * exp(d1 * tau)) / (1 - g1));
  cc1 = CC1_1 + CC1_2 * (CC1_3 - CC1_4);
  
  f1 = exp(cc1 + DD1 * v + complex(real=0,imaginary=phi * log(S)));
  
  y = Re(exp(complex(real=0,imaginary=-phi * log(K))) * f1 / (complex(real=0,imaginary=phi )));
  return(y)
}




# Trapezoidal Rule (THIS IS USED TO SLVE THE INTEGRAND)

TRAPnumint<-function(x, y)
{
  n = length(x);
  I = 0;
  for (t in 2:n) {
    I = I + 0.5 * (x[t] - x[t-1]) * (y[t-1] + y[t]);}
  return(I)
}



# Heston Option Price

Heston<-function(PutCall, kappa, theta, lambda, rho, sigma, tau, K, S, r, v)
{
  P1_int = rep(0,1001);
  P2_int = P1_int;
  phi_int = seq(0.0001,100.0001,by=.1)
  
  cnt = 1;
  
  for (phi in seq(0.0001,100.0001,by=.1))
  {
    P1_int[cnt] = HestonP1(phi, kappa, theta, lambda, rho, sigma, tau, K, S, r, v);
    P2_int[cnt] = HestonP2(phi, kappa, theta, lambda, rho, sigma, tau, K, S, r, v);
    cnt = cnt + 1;
  }
  
  p1 = 0.5 + (1 / pi) * TRAPnumint(phi_int, P1_int);
  p2 = 0.5 + (1 / pi) * TRAPnumint(phi_int, P2_int);
  
  if (p1 < 0) {p1 = 0;}
  if (p1 > 1) {p1 = 1;}
  if (p2 < 0) {p2 = 0;}
  if (p2 > 1) {p2 = 1;}
  
  
  HestonC = S * p1 - K * exp(-r * tau) * p2;
  
  if (PutCall=='Call') {y = HestonC;} else {
    if(PutCall=='Put') {y = HestonC + K * exp(-r * tau) - S}} 
  return(y)
}




#### TRYING AN EXAMPLE


K = 100;
r = .05;
rho = -.7;
kappa = 2; # lambda in our notes and in Gatheral
theta = 0.01; # v bar in our notes
sigma = 0.1; # eta in our notes
v = 0.01; # initial variance 
im_k = .5; # technical
lambda= 0.05
tau= 0.08 #time to maturity


Call_Prices = rep(1);{
Call_Prices[1]<- Heston(S=100,K=100,r=0.05,v=0.01,theta=.01,kappa=2,sigma=.1,rho=-.7, lambda=0.05, tau=.08, PutCall= 'Call');
}
Call_Prices[1]

## can either specify the variable separately, or in the function- i for some reason did both
## can extend call price by changing value in rep()

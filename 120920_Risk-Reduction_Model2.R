# Title: Consumer Risk-Reduction Behavior and New Brand Purchase - Model 2 Rev
# Date: September 20, 2019
# Econometric Model: Multiple Discrete Continuous Extreme Value Model (Bhat 2008)
# Estimation Method: Maximum Likelihood Method

# Set home directory.
setwd("C:/Users//R_MDCEV")

# Load data.
d <- read.csv('MDCEV_Input4.csv', header = T)

# Define variables
# Yoplait
p1 <- d$PR_YO/d$PR_Out
q1 <- d$OZ_YO/100
i1 <- d$ID_YO
x1 <- d$E12_YO
  
# Private brand
p2 <- d$PR_CTL/d$PR_Out
q2 <- d$OZ_CTL/100
i2 <- d$ID_CTL
x2 <- d$E12_CTL
  
# Dannon
p3 <- d$PR_DAN/d$PR_Out
q3 <- d$OZ_DAN/100
i3 <- d$ID_DAN
x3 <- d$E12_DAN
  
# Yoplait Whips!
p4 <- d$PR_YOWH/d$PR_Out
q4 <- d$OZ_YOWH/100
i4 <- d$ID_YOWH
x4 <- d$E12_YOWH

# Yoplait Light Thick & Creamy
p5 <- d$PR_YOLI/d$PR_Out
q5 <- d$OZ_YOLI/100
i5 <- d$ID_YOLI
x5 <- d$E12_YOLI

# Outside option
p0 <- d$PR_Out/d$PR_Out
q0 <- d$QT_Out/100
i0 <- d$ID_Out
  
# Demographics
si <- d$HH_Size/10
inc <- d$HH_Inco/1000
ag <- d$HH_Age/100
ep <- d$HH_Emp/100
ed <- d$HH_Ed
ed1 <- d$HH_Ed1
ed2 <- d$HH_Ed2
ed3 <- d$HH_Ed3
ed4 <- d$HH_Ed4
ed5 <- d$HH_Ed5
ed6 <- d$HH_Ed6
ed7 <- d$HH_Ed5+d$HH_Ed6
cs <- d$CONST
tr <- d$TRIPNUM
ct <- d$CATTRIP
pc <- d$Pcycle/100

# Factorial
m <- d$ID_Out+d$ID_YO+d$ID_CTL+d$ID_DAN+d$ID_YOWH+d$ID_YOLI
mf <- factorial(m-1)

# Define log likelihood function
mdcev.llf <- function(theta) {
  
  BB0 <- theta[1]
  B1 <- theta[2]
  B2 <- theta[3]
  B3 <- theta[4]
  B4 <- theta[5]
  B5 <- theta[6]
  BB1 <- theta[7]
  BB2 <- theta[8]
  BB3 <- theta[9]
  BB4 <- theta[10]
  BB5 <- theta[11]
  DD1 <- theta[12]
  DD2 <- theta[13]
  DD3 <- theta[14]
  DD4 <- theta[15]
  DD5 <- theta[16]
  KK1 <- theta[17]
  KK2 <- theta[18]
  KK3 <- theta[19]
  KK4 <- theta[20]
  KK5 <- theta[21]
  LL1 <- theta[22]
  LL2 <- theta[23]
  LL3 <- theta[24]
  LL4 <- theta[25]
  LL5 <- theta[26]
  CC1 <- theta[27]
  CC2 <- theta[28]
  CC3 <- theta[29]
  CC4 <- theta[30]
  CC5 <- theta[31]
  EE1 <- theta[32]
  EE2 <- theta[33]
  EE3 <- theta[34]
  EE4 <- theta[35]
  EE5 <- theta[36]
  
  a0 <- (1-exp(-(BB0)))
  a1 <- (1-exp(-(BB1+DD1*x1+KK1*si+LL1*inc+CC1*pc+EE1*ed7)))
  a2 <- (1-exp(-(BB2+DD2*x2+KK2*si+LL2*inc+CC2*pc+EE2*ed7)))
  a3 <- (1-exp(-(BB3+DD3*x3+KK3*si+LL3*inc+CC3*pc+EE3*ed7)))
  a4 <- (1-exp(-(BB4+DD4*x4+KK4*si+LL4*inc+CC4*pc+EE4*ed7)))
  a5 <- (1-exp(-(BB5+DD5*x5+KK5*si+LL5*inc+CC5*pc+EE5*ed7)))
  
  f0 <- (1-a0)/(q0+0)
  f1 <- (1-a1)/(q1+1)
  f2 <- (1-a2)/(q2+1)
  f3 <- (1-a3)/(q3+1)
  f4 <- (1-a4)/(q4+1)
  f5 <- (1-a5)/(q5+1)
  
  g <- ((f0^i0)*(f1^i1)*(f2^i2)*(f3^i3)*(f4^i4)*(f5^i5))
  
  h <- ((p0/f0)*i0+(p1/f1)*i1+(p2/f2)*i2+(p3/f3)*i3+(p4/f4)*i4+(p5/f5)*i5)
  
  v0 <- exp(((a0-1)*log(q0))/1)
  v1 <- exp(((B1)*i1+(a1-1)*log(q1/1+1)-log(p1))/1)
  v2 <- exp(((B2)*i2+(a2-1)*log(q2/1+1)-log(p2))/1)
  v3 <- exp(((B3)*i3+(a3-1)*log(q3/1+1)-log(p3))/1)
  v4 <- exp(((B4)*i4+(a4-1)*log(q4/1+1)-log(p4))/1)
  v5 <- exp(((B5)*i5+(a5-1)*log(q5/1+1)-log(p5))/1)
  
  vv <- (v0^i0)*(v1^i1)*(v2^i2)*(v3^i3)*(v4^i4)*(v5^i5)
  
  vt <- (v0+v1+v2+v3+v4+v5)^m 

  llf <- sum(log(sqrt((g)^2)*sqrt((h)^2)*sqrt((vv/vt)^2)*mf))
  
  return(llf)
}

# Set starting values.
theta0 <- rnorm(36,mean=0,sd=0.1)

# Maximize log likelihood function.
t0 <- proc.time()
res <- optim(theta0,mdcev.llf,gr=NULL,method="BFGS",
             hessian=TRUE,control=list(fnscale=-1,maxit=100,reltol=1e-8))
t1 <- proc.time()-t0

res$convergence
res$message
res$counts

# Obtain coefficient estimates.
theta_hat <- res$par

# Calculate standard errors.
stderr <- sqrt(-diag(solve(res$hessian)))

# Calculate t-value
tval <- theta_hat/stderr

# Obtain log likelihood value
llfv <- res$value

# Calculate AIC
aic <- -2*llfv+2*length(res$par)

# Calculate BIC
bic <- -2*llfv+log(nrow(d))*length(res$par)

# Export the Estimation Results
write.table(stderr, "C:/Users/R_MDCEV/120919_Risk-Reduction_Model2_Rev_Output.txt", 
            sep = ",", quote=F, col.names=T, row.names=F, append=F)

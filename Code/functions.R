#################################################################################
# These functions are written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
#################################################################################

#################################################################################
# error computation
#################################################################################
mbe <- function(pred, obs) {mean(pred-obs)/mean(obs)*100}
mae <- function(pred, obs) {mean(abs(pred-obs))/mean(obs)*100}
rmse <- function(pred, obs) {sqrt(mean((pred-obs)^2))/mean(obs)*100}

#################################################################################
# solar positioning
#################################################################################
calZen <- function(Tm, lat, lon, tz = 0, alt = 0)
{
  require(insol)
  
  jd <- JD(Tm)
  sunv <- sunvector(jd, lat, lon, tz)
  azi <- round(sunpos(sunv)[,1],3)#azimuth of the sun
  zen <- round(sunpos(sunv)[,2],3)#zenith angle
  #surface.norm = normalvector(tilt, orientation)
  #inc = round(as.numeric(degrees(acos(sunv%*% as.vector(surface.norm)))),3)
  doy <- daydoy(Tm)
  da <- (2 * pi / 365) * (doy - 1)
  re = 1.000110+0.034221*cos(da)+0.001280*sin(da)+0.00719*cos(2*da)+0.000077*sin(2*da)
  E0n = round(1361.1*re,3)#extraterrestrial direct normal irradiance
  E0 = round(1361.1*re*cos(radians(zen)))#horizontal extraterrestrial irradiance
  E0 <- ifelse(zen>=90, 0, E0)
  
  out = list(zen, azi, E0)
  names(out) = c("zenith", "azimuth", "E0")
  out
}


#################################################################################
# separation models
#################################################################################
Erbs <- function(Kt)
{
  if(Kt <= 0.22)
  {
    Kd <- 1- 0.09*Kt
  }else if(Kt > 0.22 & Kt <= 0.80){
    Kd <- 0.9511 - 0.1604*Kt + 4.388*(Kt)^2 - 16.638*(Kt)^3 + 12.336*(Kt)^4
  }else{
    Kd <- 0.165
  }
  Kd
}

Engerer2 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  Z <- data$Z
  delta_ktc <- data$delta_ktc
  kde <- data$kde
  # model
  C <- 4.2336e-2; beta0 <- -3.7912; beta1 <- 7.5479; beta2 <- -1.0036e-2; beta3 <- 3.1480e-3; beta4 <- -5.3146; beta5 <- 1.7073;
  C + (1-C)/(1+exp(beta0+beta1*kt + beta2*AST + beta3*Z + beta4*delta_ktc)) + beta5*kde
}

Engerer4 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  Z <- data$Z
  delta_ktc <- data$delta_ktc
  kde <- data$kde
  # model
  C <- 1.0562e-1; beta0 <- -4.1332; beta1 <- 8.2578; beta2 <- 1.0087e-2; beta3 <- 8.8801e-4; beta4 <- -4.9302; beta5 <- 4.4378e-1;
  C + (1-C)/(1+exp(beta0+beta1*kt + beta2*AST + beta3*Z + beta4*delta_ktc)) + beta5*kde
}

Starke1 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  psi <- data$psi
  Gcsky <- data$Gcsky/277.78  # coefficient fitted using MJ/h/m2
  kcsi <- data$kcsi
  # model
  beta <- c(-6.70407, 6.99137, -0.00048, 0.03839, 3.36003, 1.97891, -0.96758, 0.15623, -4.21938, -0.00207, -0.06604, 2.12613, 2.56515, 1.62075) #australia
  #beta <- c(-6.37505, 6.68399, 0.01667, 0.02552, 3.32837, 1.97935,-0.74116, 0.19486, -3.52376, -0.00325, -0.03737, 2.68761, 1.60666, 1.07129) #brazil
  ifelse((kcsi >=1.05 & kt>0.65), 1/(1+exp(beta[8] + beta[9]*kt + beta[10]*AST + beta[11]*alpha + beta[12]*daily_kt + beta[13]*psi + beta[14]*Gcsky)), 1/(1+exp(beta[1] + beta[2]*kt + beta[3]*AST + beta[4]*alpha + beta[5]*daily_kt + beta[6]*psi + beta[7]*Gcsky)) )
}

Starke2 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  psi <- data$psi
  Gcsky <- data$Gcsky/277.78  # coefficient fitted using MJ/h/m2
  kcsi <- data$kcsi
  # model
  #beta <- c(-6.70407, 6.99137, -0.00048, 0.03839, 3.36003, 1.97891, -0.96758, 0.15623, -4.21938, -0.00207, -0.06604, 2.12613, 2.56515, 1.62075) #australia
  beta <- c(-6.37505, 6.68399, 0.01667, 0.02552, 3.32837, 1.97935,-0.74116, 0.19486, -3.52376, -0.00325, -0.03737, 2.68761, 1.60666, 1.07129) #brazil
  ifelse((kcsi >=1.05 & kt>0.65), 1/(1+exp(beta[8] + beta[9]*kt + beta[10]*AST + beta[11]*alpha + beta[12]*daily_kt + beta[13]*psi + beta[14]*Gcsky)), 1/(1+exp(beta[1] + beta[2]*kt + beta[3]*AST + beta[4]*alpha + beta[5]*daily_kt + beta[6]*psi + beta[7]*Gcsky)) )
}

Starke3 <- function(data, clim)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  hourly_kt <- data$hourly_kt
  psi <- data$psi
  Gcsky <- data$Gcsky#/277.78  # coefficients seem to be fitted using W/m2
  kcsi <- data$kcsi
  clim <- substr(clim, 1, 1)
  clim <- ifelse(clim %in% c("A", "B", "C", "D", "E"), clim, "C")
  # model
  clims <- c("A", "B", "C", "D", "E")
  beta0 <- switch(match(clim, clims), 0.29566, -1.7463, -0.083, 0.67867, 0.51643)
  beta1 <- switch(match(clim, clims), -3.64571, -2.20055, -3.14711, -3.79515, -5.32887)
  beta2 <- switch(match(clim, clims), -0.00353, 0.01182, 0.00176, -0.00176, -0.00196)
  beta3 <- switch(match(clim, clims), -0.01721, -0.03489, -0.03354, -0.03487, -0.07346)
  beta4 <- switch(match(clim, clims), 1.7119, 2.46116, 1.40264, 1.33611, 1.6064)
  beta5 <- switch(match(clim, clims), 0.79448, 0.70287, 0.81353, 0.76322, 0.74681)
  beta6 <- switch(match(clim, clims), 0.00271, 0.00329, 0.00343, 0.00353, 0.00543)
  beta7 <- switch(match(clim, clims), 1.38097, 2.30316, 1.95109, 1.82346, 3.53205)
  beta8 <- switch(match(clim, clims), -7.00586, -6.53133, -7.28853, -7.90856, -11.70755)
  beta9 <- switch(match(clim, clims), 6.35348, 6.63995, 7.15225, 7.63779, 10.8476)
  beta10 <- switch(match(clim, clims), -0.00087, 0.01318, 0.00384, 0.00145, 0.00759)
  beta11 <- switch(match(clim, clims), 0.00308, -0.01043, 0.02535, 0.10784, 0.53397)
  beta12 <- switch(match(clim, clims), 2.89595, 1.73562, 2.35926, 2.00908, 1.76082)
  beta13 <- switch(match(clim, clims), 1.13655, 0.85521, 0.83439, 1.12723, 0.41495)
  beta14 <- switch(match(clim, clims), -0.0013, -0.0003, -0.00327, -0.00889, -0.03513)
  beta15 <- switch(match(clim, clims), 2.75815, 2.63141, 3.19723, 3.72947, 6.04835)
  ifelse((kcsi >=1.05 & kt>0.75), 1/(1+exp(beta0 + beta1*kt + beta2*AST + beta3*alpha + beta4*daily_kt + beta5*psi + beta6*Gcsky + beta7*hourly_kt)), 1/(1+exp(beta8 + beta9*kt + beta10*AST + beta11*alpha + beta12*daily_kt + beta13*psi + beta14*Gcsky + beta15*hourly_kt)) )
}

Abreu <- function(data, clim)
{
  # get input from tibble
  kt <- data$kt
  clim <- substr(clim, 1, 1)
  if(clim == "A")
  {
    clim <- "TR"
  }else if(clim == "B"){
    clim <- "AR"
  }else if(clim == "E"){
    clim <- "HA"
  }else{
    clim <- "TM"
  }
  # model
  clims <- c("AR", "HA", "TM", "TR")
  A <- switch(match(clim, clims), 11.39, 7.83, 10.79, 11.59)
  B <- switch(match(clim, clims), -6.25, -4.59, -5.87, -6.14)
  n <- switch(match(clim, clims), 1.86, 3.25, 2.24, 1.87)
  (1+(A*(kt-0.5)^2+B*(kt-0.5)+1)^(-n))^(-1/n)
}

Paulescu <- function(data)
{
  # get input from tibble
  kt <- data$kt
  daily_kt <- data$daily_kt
  # model
  1.0119 - 0.0316*kt - 0.0294*daily_kt -1.6567*(kt-0.367)*ifelse(kt-0.367>=0, 1, 0) + 1.8982*(kt-0.734)*ifelse(kt-0.734>=0, 1, 0) -0.8548*(daily_kt-0.462)*ifelse(daily_kt-0.462>=0, 1, 0)
}

Every1 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  psi <- data$psi
  # model
  beta0 <- -6.862; beta1 <- 9.068; beta2 <- 1.468e-2; beta3 <- -4.72e-3; beta4 <- 1.703; beta5 <- 1.084;
  1/(1+exp(beta0+beta1*kt + beta2*AST + beta3*alpha + beta4*daily_kt + beta5*psi)) 
}

Every2 <- function(data, clim)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  psi <- data$psi
  clim <- ifelse(clim %in% c("Am", "Aw", "BSh", "BSk", "BWh", "Cfa", "Cfb", "Csa", "Csb"), clim, "Other")
  # model
  # some climates are not included by Every2, so the orignal BRL is used instead
  clims <- c("Other", "Am", "Aw", "BSh", "BSk", "BWh", "Cfa", "Cfb", "Csa", "Csb")
  beta0 <- switch(match(clim, clims), -5.38, -6.433, -6.047, -6.734, -7.310, -7.097, -6.484, -6.764, -7.099, -7.080)
  beta1 <- switch(match(clim, clims), 6.63, 8.774, 7.540, 8.853, 10.089, 9.416, 8.301, 9.958, 10.152, 10.460)
  beta2 <- switch(match(clim, clims), 6e-3, -4.4e-4, 6.24e-3, 2.454e-2, 1.852e-2, 1.254e-2, 1.577e-2, 1.271e-2, -2.6e-4, 9.64e-3)
  beta3 <- switch(match(clim, clims), -7e-3, -5.78e-3, -2.99e-3, -4.95e-3, -6.93e-3, -4.16e-3, -3.38e-3, -1.249e-2, -7.44e-3, -1.420e-2)
  beta4 <- switch(match(clim, clims), 1.75, 2.096, 2.077, 1.874, 1.296, 1.661, 1.607, 0.928, 1.147, 1.134)
  beta5 <- switch(match(clim, clims), 1.31, 0.684, 1.208, 0.939, 1.114, 1.130, 1.307, 1.142, 1.184, 1.017)
  1/(1+exp(beta0+beta1*kt + beta2*AST + beta3*alpha + beta4*daily_kt + beta5*psi)) 
}

Yang4 <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  Z <- data$Z
  delta_ktc <- data$delta_ktc
  kde <- data$kde
  kds <- data$kds60 # 1-h diffuse fraction from Engerer2
  clim <- ifelse(clim %in% c(1,2,3,4,5,6), clim, "Other")
  # model
  C <- 0.0361; beta0 <- -0.5744; beta1 <- 4.3184; beta2 <- -0.0011; beta3 <- 0.0004; beta4 <- -4.7952; beta5 <- 1.4414; beta6 <- -2.8396;
  C + (1-C)/(1+exp(beta0+beta1*kt + beta2*AST + beta3*Z + beta4*delta_ktc+ beta6*kds)) + beta5*kde 
}


Yang5 <- function(data, clim) # same as Yang4 but regime-switching
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  Z <- data$Z
  delta_ktc <- data$delta_ktc
  kde <- data$kde
  kds <- data$kds60 
  # model
  clims <- c(1,2,3,4,5,6)
  C <- switch(match(clim, clims), 0.131047425,-0.01614172,-2.747523e-01,-0.010948960,0.042971966)
  beta0 <- switch(match(clim, clims), -4.267395518,-3.33037519,3.608535e-01,-0.921287782,-1.644372959)
  beta1 <- switch(match(clim, clims), 7.680511140,5.72307441,3.986924e-01,3.650149139,4.718078475)
  beta2 <- switch(match(clim,clims), 0.005399573,0.01295824,4.789620e-03,0.007674206,0.014623824)
  beta3 <- switch(match(clim,clims), 0.017484478,0.01230061,3.871555e-04,0.004936045,0.007453108)
  beta4 <- switch(match(clim,clims), 0.915902661,-0.96482793,-1.020264e+01,-3.764652421,-3.352233222)
  beta5 <- switch(match(clim,clims), 0.521761178,0.94203658,2.124752e+00,1.364819177,1.251921688)
  beta6 <-switch(match(clim,clims), -1.688185184,-1.68332469,-1.784547e+00,-2.118672147,-2.364771589)
  
  C + (1-C)/(1+exp(beta0+beta1*kt + beta2*AST + beta3*Z + beta4*delta_ktc+ beta6*kds)) + beta5*kde 
}

BRL <- function(data)
{
  # get input from tibble
  kt <- data$kt
  AST <- data$AST
  alpha <- 90-data$Z
  daily_kt <- data$daily_kt
  psi <- data$psi
  # model
  beta <- c(-5.38, 6.63, 0.006, -0.007, 1.75, 1.31)
  1/(1+exp(beta[1] + beta[2]*kt + beta[3]*AST + beta[4]*alpha + beta[5]*daily_kt + beta[6]*psi))
}




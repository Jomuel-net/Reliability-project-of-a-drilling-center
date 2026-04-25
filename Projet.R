# imports
library (survival)
library(flexsurv)
library(fitdistrplus)
library(phasetype)

# ## ----------------- question 2 ( estimation de la loi de T de la plateforme rotative) -----------------

df.temps.avant.defaillance <- read.csv("Temps_avant_defaillance.csv", sep = ';')
head(df.temps.avant.defaillance)

# objet de survi

objet.survi <- Surv (time = df.temps.avant.defaillance$time, event = df.temps.avant.defaillance$status)

# kaplan Maier (fonction de survi du dispositif) [ sans connaitre la loi]

km <- survfit (objet.survi~1)
plot (km, xlab = "time", ylab = "survie", main ="kaplan Maier" ) # fournit plus une courbe concave : donc weibull avec beta < 1?


# estimation

# test des lois possibles
test.expo <- flexsurvreg(objet.survi~1, dist = "exponential")
test.weibull <- flexsurvreg(objet.survi~1, dist = "weibull")
test.gamma <- flexsurvreg(objet.survi~1, dist = "gamma")
test.lognormal <- flexsurvreg(objet.survi~1, dist = "lognormal")


print(test.expo)
print(test.weibull)   # weibull est correcte
print(test.gamma)
print(test.lognormal)

lines( test.weibull , col = "red")

#----------------------------------------Question 3 (loi du sys de lubrification)----------------------------------------------------

df.TTF.lubr <- read.csv("TTF_sys_lubr.csv", 
                        sep = ";", 
                        dec = ",")

head(df.TTF.lubr)


# estimation ( sans censure)
ttf <- df.TTF.lubr$TTF
fit_weib <- fitdist(ttf, "weibull")
print(fit_weib)  

#------------------------------------Question 4 ( loi du temps complet et partiel)---------------------------------------------- 

                                      #-----temps complet-----
df.Tcomplet.lubr <- read.csv("Tcomplet.csv", sep = ';', dec = ",")
head(df.Tcomplet.lubr)

# estimation ( sans censure)

tcomplet <- df.Tcomplet.lubr$Tcomplet

# virgule par un point 

tcomplet <- gsub(",", ".", tcomplet)

# convertir en numérique
tcomplet <- as.numeric(tcomplet)

fit_lognorm <- fitdist(tcomplet, distr = "lnorm")
print(fit_lognorm)  


                                      #------temps partiel-------
df.Tpartiel.lubr <- read.csv("Tpartiel.csv", sep = ';', dec = ",")
head(df.Tpartiel.lubr)


# estimation ( sans censure)
tpartiel <- df.Tpartiel.lubr$Tpartiel


# virgule par un point 

tpartiel <- gsub(",", ".", tpartiel)

# convertir en numérique
tpartiel <- as.numeric(tpartiel)

fit_lognorm <- fitdist(tpartiel, distr="lnorm")
print(fit_lognorm)  



#------------------------------------Question 5 (dispo moyenne du sys)---------------------------------------------- 
# dispo moyenne = MTTF / MTBF

ttf.ttr.df <- read.csv("TTF+TTR_sys_lubr.csv", sep = ';')
MTTF <- mean(df.TTF.lubr$TTF)
MTBF <- mean (ttf.ttr.df$TTF.TTRpar)
print (MTTF/MTBF)

#------------------------------------Question 7 (impact sys lubr sur meche)---------------------------------------------- 
# M = A*Q + (1-A)Q'



#-----------------------------Maintenance-------------------------------------------------------------
#------------------------systeme de lubrification-------------
# maintenance actuelle (ok)
# maintenance périodique
# estimation de la convolution de la fonction de repartition du sys de lubr (monte-carlo)


M_function <- function(b) {
  
  alpha <- 5819.75      # scale
  beta <- 2.23          # shape
  kmax <- 20            # ordre max de convolution
  N <-  1e5             # taille de simulation pour chaque k
  
  somme_Fk <- 0
  for (k in 1:kmax) {
    S_k <- replicate(N, sum(rweibull(k, shape = beta, scale = alpha)))
    Fk_hat <- mean(S_k <= b)   # estimation de F^{*k}(b)
    somme_Fk <- somme_Fk + Fk_hat
  }
  
  return(somme_Fk)        # approx de sum_{k=1}^∞ F^{*k}(b)
  
}

# cc (MTTR du sys lubr) et cp
MTTR <- MTBF - MTTF
cc <- MTTR * 100 +250 
cp <- 4*100 + 900 

# optimisation du coût par unité de temps
cout_optim <- function(b) {
  return ((cc * M_function(b) +cp)/b)
}

# generation de b_0* (h)
optim <- optimize(cout_optim, interval = c (300,9000), tol = 1e-8 )  # interval autour du scale 
print(optim$minimum)     # 8996.73h
print(optim$objective)   # 0.2 et donc 0.2*8996.73 = 1300 sur un cyle 


#-------------------------------Bloc motorisation------------------------
# maintenance actuelle (ok)
# maintenance périodique

M_function <- function(b) {
  
  alpha <- 13439.77     # scale
  beta <- 2.35          # shape
  kmax <- 20            # ordre max de convolution
  N <-  1e5             # taille de simulation pour chaque k
  
  somme_Fk <- 0
  for (k in 1:kmax) {
    S_k <- replicate(N, sum(rweibull(k, shape = beta, scale = alpha)))
    Fk_hat <- mean(S_k <= b)   # estimation de F^{*k}(b)
    somme_Fk <- somme_Fk + Fk_hat
  }
  
  return(somme_Fk)        # approx de sum_{k=1}^∞ F^{*k}(b)
  
}

# MTTR,  cc et cp
MTTR <- 2.5
cc <- 1600
cp <- 1400

# optimisation du coût par unité de temps
cout_optim <- function(b) {
  return ((cc * M_function(b) +cp)/b)
}

# generation de b_0* (h)
optim <- optimize(cout_optim, interval = c (10000,20000), tol = 1e-8 )  # interval autour du scale 
print(optim$minimum)     # 19908.33h
print(optim$objective)   # 0.17 et donc 3384,4161  sur un cyle 

#----------------------- Maintenance par âge motorisation -----------------------

# fonction tau
tau <- function(a0) {
  alpha <- 13439.77     # scale
  beta <- 2.35          # shape
  # Intégrale de 0 à a0 de u f(u) du
  I <- integrate(function(u) u * dweibull(u, shape = beta, scale = alpha),lower = 0, upper = a0)$value
  
  # Terme a0 * R(a0)
  T <- a0 * pweibull(a0, shape = beta, scale = alpha, lower.tail = FALSE)
  
  return(I + T)
}

# optimisation du coût par unité de temps
C_A <- function(a0) {
  alpha <- 13439.77     # scale
  beta <- 2.35          # shape
  (cc * pweibull(a0, shape = beta, scale = alpha) +
     cp * pweibull(a0, shape = beta, scale = alpha, lower.tail = FALSE)) /
    tau(a0)
}

# génération de a0* (h)
optim <- optimize(C_A, interval = c(10000, 15000), tol = 1e-8)
print(optim$minimum)
print(optim$objective)





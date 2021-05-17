rm(list=ls())
require(zoo) # Pour formalise les "time series" de manières pratiques
require(tseries) # Pour diverses fonctions liées au "time series"
require(fUnitRoots) # Pour le test de racine unitaire avec adfTest
require(latex2exp) # Pour ecrire du texte en LaTeX
require(forecast) # Pour la prevision


###########
# Partie 1#
###########

# On change le repertoire de travail
#setwd(dir="C:/Users/roisc/OneDrive/Documents/Mes projets/Projet Séries Temporelles linéaires")
# On importe la série la base de données au format .csv
mydatabase <- read.csv(file = file.choose(), header = TRUE, sep = ";", dec = ",")
# On formalise notre TS (TS pout "time serie")
date_char <- as.character(mydatabase$date)
date_char[1];date_char[length(date_char)]
date.int <- as.yearmon(seq(from=1990+0/12,to=2021+2/12,by=1/12))
xt.init <- zoo(mydatabase[,2], order.by=date.int)
# On met dans T.int la longueur de la TS initial
T.int <- length(xt.init)
# On enleve les deux dernieres valeurs, pour pouvoir comaparer avec les predictions
xt <- xt.init[1 : (T.int-2)] 
date <- date.int[1 : (T.int-2)]
# On met dans T la longueur de la TS qu'on utilisara dans la suite
T <- length(xt)
# On plot la TS xt
plot(xt, xlab = "dates", ylab = TeX("$x_t$"), main = TeX("$x_t =$ Indice de production industrielle"))
# La variabilite est un peu trop grande. On va la reduire
# On va donc opter pour une linearisation de la serie, par transformation logarithmique
yt <- log(xt)
plot(yt, xlab = "dates", ylab = TeX("$y_t$"), main = TeX("$y_t = log(x_t)$"))
# Elle semble a present ne plus avoir une trop grande variabilite. 
# Neanmoins, elle semble osciller autour d'une constante.Regardons cela de plus pres.
summary(lm(yt~1))
# On voit que la constante est tres significative. Nos soupçons semblent donc averes.
# On va faire un teste de racine unitaire pour voir si elle est stationnaire autour de cette constante.
# Cette fonction permettant de realiser le test d'autocorrelation de Ljung-Box sur une TS
Qtests <- function(series, k, fitdf=0){
  pvals <- apply(matrix(1:k), 1, FUN=function(l){
    pval <- if (l<=fitdf) NA else Box.test(x = series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
# Cette fonction realise le test de DF ou ADF a partir du moment ou les residus de la regression 
# presentent une abscence d'autocorrelation
adfTest_valid <- function(series, kmax, type){
  k <- 0
  noautocorr <- 0
  while(noautocorr==0 & k<=kmax){
    cat(paste0("ADF with ",k, " lags : residuals OK ? "))
    adf <- adfTest(series, lags=k, type=type)
    pvals <- Qtests(adf@test$lm$residuals,kmax,length(adf@test$lm$coefficients))[,2]
    if(sum(pvals<0.05,na.rm=T) == 0){
      noautocorr <- 1; cat("YES \n")
    }
    else cat("NO \n")
    k <- k + 1
  }
  return(adf)
}
# On realise donc le test de DF ou ADF (selon que le lag vaut 0 ou pas)
adfTest_valid(yt, 36, "c")
# Le lag a partir duquel les residus ne sont plus autocorreles est 4.
# La p_value est < 0,1. Donc on rejette l'hypothese de presence de racine unitaire a 1%. 
#La serie yt est donc stationnaire autour de sa moyenne. On va retirer cette moyenne et retester.
zt = yt - mean(yt)
adfTest_valid(zt, 36, "nc")
# zt est donc stationnaire, la p_value est < 0,1.

###########
# Partie 2#
###########

# ACF et PACF
nlag = 3*12
acf=Acf(zt,lag.max=nlag,plot=FALSE)
pacf=Pacf(zt,lag.max=nlag,plot=FALSE)
autoplot(acf, lwd=2, xlab="Retards", ylab="ACF", main=TeX("ACF of $z_t$"))
autoplot(pacf, lwd=2, xlab="Retards", ylab="ACF", main=TeX("ACF of $z_t$"))

# On voit que pmax = 6 et qmax = 18

# Fonction de test des significativites individuelles des coefficients
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=zt, k=24){
estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
checks <- c(arsignif,masignif,resnocorr)
ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}
## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax = 3 et q<=max = 73
armamodelchoice <- function(pmax,qmax){
pqs <- expand.grid(0:pmax,0:qmax) 
m <- t(apply(matrix(1:dim(pqs)[1]),1,function(row){
                                          p <- pqs[row,1] 
                                          q <- pqs[row,2]
                                          cat(paste0("Computing ARMA(",p,",",q,") \n"))
                                          modelchoice(p,q)
                                          }))
return(m)
}
armamodels <- armamodelchoice(6,18) #estime tous les arima (patienter...)
selec <- armamodels[armamodels[,"ok"]==1 & !is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
#creation d'une liste des ordres p et q des modeles candidats
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les elements de la liste
models <- lapply(pqs, function(pq) arima(zt,c(pq[["p"]],0,pq[["q"]]))) #creation une liste des modeles
#calcule les AIC et BIC des modeles candidats
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m)))
# On va utiliser le BIC pour pouvoir plus penaliser les modeles trop complexes.
# On voit que par BIC le modele retenu est le ARMA(1,2). On l'estime.
## Fonction d'affichage des tests pour la sélection du modèle ARIMA
arimafit <- function(estim){
  adjust <- round(signif(estim), 5)
  pvals <- Qtests(estim$residuals, 24, length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat(" tests de nullité des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrélation des résidus : \n")
  print(pvals)
}
# ARIMA(1,0,2)
arima102 <- arima(zt,c(1,0,2)) # Estimation du ARMA(1,2)
#signif(arima102)
arimafit(estim = arima102)
# L'intercept est significativement nul, ce qui est naturel puisque zt est centre.
# Les trois autres coefficients sont tres significatives, le modele est bien ajuste.
# Le modele passe les tests portemanteaux, il est donc bien valide.

# Normalite des residus
jarque.bera.test(arima102$residuals)
qqnorm(arima102$residuals,main="QQ-plot des rÃ©sidus")
qqline(arima102$residuals)
# Heteroscedasticite des residus
t = 1:T
bptest(arima102$residuals ~ t)


###########
# Partie 3#
###########

# Prevision

pred_arima102 <- forecast(arima102, h = 2, level = 95)
autoplot(pred_arima102, include = 30, PI = TRUE, lwd = 2)
pred_xt_m = exp(mean(yt) + pred_arima102$mean)
pred_xt_l = exp(mean(yt) + pred_arima102$lower)
pred_xt_u = exp(mean(yt) + pred_arima102$upper)
tab <- cbind(pred_xt_m, pred_xt_l, pred_xt_u) 
tab


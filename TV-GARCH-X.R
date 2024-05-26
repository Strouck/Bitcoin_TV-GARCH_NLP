################################################################################
#                                 GARCH                                        #
################################################################################
################################################################################
#                                Libraries                                     #
################################################################################

library(quantmod)
library(dplyr)
library(tidyverse)
library(tseries)
library(xts)
library(readr)
library(moments)
library(forecast)
library(rugarch)
library(MCS)
library(Hmisc)
library(aTSA)
library(FinTS)
library(pracma)

################################################################################
#                                 Dataset                                      #
################################################################################
data <- read.csv("/Users/***/Downloads/Dataset_final.csv")

################################################################################
#                             Preparing data                                   #
################################################################################
#In this work i will use exogenous variables so the sutable models for me are
#ARIMAX/ARMAX model where X stands for exogenous variables. I'll start with it!

data$Date <- as.Date(data$Date)
rownames(data) <- data$Date


#clear some columns

data$index <- NULL
data$X <- NULL

#Making some basic plots

plot(data$Date,data$Bitcoin_price, main = 'Bitcoin',type="l")
plot(data$Date,data$Gold_price, main = 'Gold',type="l")

################################################################################
#                           Investigating data                                 #
################################################################################

#######################################
#           Bitcoin + Gold            #
#######################################

#Firstly i need to compute returns

#Bitcoin_ret_log <- diff(log(data$Bitcoin_price))[-1] * 100
#Gold_ret_log <- diff(log(data$Gold_price))[-1] * 100

Bitcoin_ret <- Delt(data$Bitcoin_price)[-1]*100
Gold_ret <- Delt(data$Gold_price)[-1]*100

#Bitcoin_ret_log <- diff(log(data$Bitcoin_price))[-1] * 100
#Gold_ret_log <- diff(log(data$Gold_price))[-1] * 100

data <- data.frame(data[2:3282,])

#data <- data[2:2209,]

#Sent_crypto <- movavg(data$Sentiment_VADER_comp_crypto, n=4)
#Sent_classic <- movavg(data$Sentiment_VADER_comp_classic, n=4)

#data <- data[2:2209,]

data$Bitcoin_ret <- Bitcoin_ret
data$Gold_ret <- Gold_ret
data$EPU_index <- EPU_index

plot(data$Date,data$Bitcoin_ret, main = 'Bitcoin returns',type="l")
plot(data$Date,data$Gold_ret, main = 'Gold return',type="l")
plot(data$Date,data$VIX, main = 'VIX',type="l")
plot(data$Date,data$Sentiment_VADER_comp_crypto, main = 'Crypto sentiment',type="l")
plot(data$Date,data$Sentiment_VADER_comp_classic, main = 'Classic sentiment',type="l")

#plot(data$Date,data$Crypto_set, main = 'Crypto',type="l")
#plot(data$Date,data$Classic, main = 'Classic',type="l")

data <- data[data$Date >= '2017-01-01',]

#######################################
#             Stat. tests             #
#######################################

#Checking normality

#H0: data is normally distributed
#H1: data is not came from a normal distribution

jarque.bera.test(Bitcoin_ret)

#Returns are not normally distributed

#Plot histogram of returns

hist(Bitcoin_ret, freq = F, breaks = 'FD')
curve(dnorm(x, mean = mean(Bitcoin_ret), sd = sd(Bitcoin_ret)), 
      col="blue", lwd=2, add=TRUE)
legend("topright", legend = "Normal Distribution",
       col = "blue", lty=1:2, cex = 0.8)

#Checking stationarity №1

#H0: The time series is regarded as non-stationary
#H1: The time series is regarded as stationary

adf.test(data$Bitcoin_ret)

#So the series is stationary

# Checking stationarity №2

#H0: The series is stationary, meaning that it does not have a unit root
#H1: The series is non-stationary due to a unit root

kpss.test(data$Bitcoin_ret)

#So based on p-value we can't reject the null

matrix=matrix(c(data$Sentiment_VADER_comp_crypto,
                data$Sentiment_VADER_comp_classic,
                data$Gold_ret,
                data$VIX),ncol=4)

cor(matrix)

#######################################
#             ARIMA model             #
#######################################

#just for investigation and testing

model.arima.bitcoin <- auto.arima(data$Bitcoin_ret_log ,stepwise = FALSE)

summary(model.arima.bitcoin)

acf(model.arima.bitcoin$residuals^2)

#We have some ARCH effect. Let's check it using 
Box.test(model.arima.bitcoin$residuals^2, type = "Ljung-Box")
#H_0: The data are independently distributed

arimaModel_1=arima(data$Bitcoin_ret, order=c(5,0,0))


arch.test(arimaModel_1)

################################################################################
#                           Building GARCH model                               #
################################################################################

#######################################
#           sGARCH model             #
#######################################

# ARIMA + GARCH
#Specification with crypto specific news

garch.spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                               garchOrder = c(1, 1)), 
                         mean.model = list(armaOrder = c(5,0)), 
                         distribution.model = "snorm")
garch.model <- ugarchfit(spec = garch.spec, data = data$Bitcoin_ret)
garch.model

#Results of regarch are unstable

################################################################################
#                  Building TV - GARCH - X model                               #
################################################################################
install.packages('tvgarch',version='2.4')
library(tvgarch)

data <- data[data$Date >= '2017-01-01',]

#First of all lest again check the graph for bitcoin returns for the presence of
#volatility clustering. So the unconditional variance is not constant !

plot(data$Date,data$Bitcoin_ret, main = 'Bitcoin returns',type="l")

#Before estimating the model we need to identify the number of transitions or 
#changes in unconditional variance.

transitions <- tvgarchTest(data$Bitcoin_ret)
transitions

#so we that 1 transition is appropriate. Proceed to estimation

matrix=matrix(c(data$VIX_ret,
                data$Sentiment_VADER_comp_crypto),ncol=2)

tvgarchEst <- tvgarch(y = data$Bitcoin_ret, 
                      order.g = 1,
                      order.h = c(1,1,0),
                      xreg=matrix)
summary(tvgarchEst)

#Now we make plots for TV and Garch components
plot(tvgarchEst)

#Extract the date of transition
coef.g <- coef(tvgarchEst, spec = "tv")
c11 <- coef.g["location1"]
stime2 <- tvgarchEst$xtv
stime2[which.min(abs(stime2 - c11))]

data$Date[1819]

#Now we need to calculate the p-values for our coefficient estimates
#we need estimate, standard error and number of degrees of freedom

t_stats_TV <- tvgarchEst[["par.g"]]/tvgarchEst[["se.g"]]
t_stats_GARCH <- tvgarchEst[["iter.fit.h"]][["par"]]/tvgarchEst[["se.h"]]

p_value_TV=2*pt(abs(t_stats_TV), 
           nobs(tvgarchEst)-length(coef(tvgarchEst)),lower.tail=FALSE)

p_value_GARCH=2*pt(abs(t_stats_GARCH), 
                nobs(tvgarchEst)-length(coef(tvgarchEst)),lower.tail=FALSE)

# P-value for the TV components
p_value_TV

# P-value for GARCH components
p_value_GARCH

logLik(tvgarchEst)

toLatex(tvgarchEst)

BIC2 <- BIC(logLik(tvgarchEst))
BIC2
acf(tvgarchEst$residuals^2)

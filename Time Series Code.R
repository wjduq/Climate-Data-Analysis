## Code for Time Series Analysis

#Import Libraries
require(mosaic)
options(digits=3)
palette = trellis.par.get()$superpose.symbol$col
library(Stat2Data)
library(forecast)
library(lmtest) # if you want p-values for coefficients in ARIMA models
library(tseries) # need for Augmented Dickey-Fuller test
require(fBasics) # need for normality tests of residuals

data("CO2Hawaii")
head(CO2Hawaii)
tail(CO2Hawaii)

#Create Time Series
class(CO2Hawaii)
co2 = ts(CO2Hawaii$CO2)
class(co2)

plot(co2, ylab = "CO2 levels", xlab = "Time since January 1988", main = "CO2 Levels", col = "blue")
plot(as.ts(co2[324:360]), ylab = "CO2 levels", xlab = "Time since January 2015", main = "CO2 levels over the last 3 years in Hawaii", col = "blue")

modT=lm(CO2~t,data=CO2Hawaii)
summary(modT)

plot(modT, which = 1)

modTQ=lm(CO2~t+I(t^2),data=CO2Hawaii)
summary(modTQ)
plot(modTQ, which = 1)

plot(as.ts(modT$residuals))
plot(as.ts(modTQ$residuals))

CO2Hawaii$Xcos=cos(2*pi*CO2Hawaii$t/12)
CO2Hawaii$Xsin=sin(2*pi*CO2Hawaii$t/12)

modCos=lm(CO2~Xcos+Xsin+t+I(t^2),data=CO2Hawaii)
summary(modCos)

modseason=lm(CO2~as.factor(Month)+t+I(t^2),data=CO2Hawaii)
summary(modseason)

fun1 <- makeFun(modT)
fun2 <- makeFun(modTQ)
fun3 <- makeFun(modCos)
fun4 <- makeFun(modseason)

fun1(t = 370)
fun2(t = 370)
fun3(cos(2*pi*370/12),sin(2*pi*370/12),370)
fun4(10,t = 370)

AIC(modCos)
AIC(modseason)

real2018 = c(407.98,
             408.36,
             409.21,
             410.24,
             411.23,
             410.81,
             408.83,
             407.02,
             405.52,
             405.93,
             408.04,
             409.17)

cos2018 <- c(fun3(cos(2*pi*361/12),sin(2*pi*361/12),361), fun3(cos(2*pi*362/12),sin(2*pi*362/12),362), fun3(cos(2*pi*363/12),sin(2*pi*363/12),363), fun3(cos(2*pi*364/12),sin(2*pi*364/12),364), fun3(cos(2*pi*365/12),sin(2*pi*365/12),365), fun3(cos(2*pi*366/12),sin(2*pi*366/12),366), fun3(cos(2*pi*367/12),sin(2*pi*367/12),367), fun3(cos(2*pi*368/12),sin(2*pi*368/12),368), fun3(cos(2*pi*369/12),sin(2*pi*369/12),369), fun3(cos(2*pi*370/12),sin(2*pi*370/12),370), fun3(cos(2*pi*371/12),sin(2*pi*371/12),371), fun3(cos(2*pi*372/12),sin(2*pi*372/12),372))

seasons2018 <- fun4(c(1:12), t = c(361:372))

(cor(real2018, cos2018))^2
(cor(real2018, seasons2018))^2

modseason2=lm(CO2~as.factor(Month),data=CO2Hawaii)
summary(modseason2)

season2fit <- ts(data.frame(co2, modseason$fitted))
plot(season2fit,plot.type="s",lwd=3,col=c("red","purple"),ylab="CO2 levels",main="CO2 levels with seasonal fit (modseason)")

season2fit <- ts(data.frame(co2, modseason2$fitted))
plot(season2fit,plot.type="s",lwd=3,col=c("red","purple"),ylab="CO2 levels",main="CO2 levels with seasonal fit (modseason2)")

par(mfrow=c(2,2))
acf(modseason$residuals)
pacf(modseason$residuals)

mod200 = Arima(modseason$residuals,order=c(2,0,0))

acf(mod200$residuals)
pacf(mod200$residuals)
Box.test(mod200$residuals, type = "Ljung-Box")
adf.test(mod200$residuals)

coeftest(mod200)

acf(mod200$residuals)
Box.test(mod200$residuals, type = "Ljung-Box")

tsdiag(mod200)

adf.test(modseason$residuals)
Box.test(modseason$residuals, type = "Ljung-Box")
adf.test(mod200$residuals)

pred5=forecast(mod200,h=12)
fun1 <- makeFun(modseason)
ft0 = fun1(t = c(361:372), Month = c(1:12))
pred5list <- c(0.0638,0.0681, 0.0577, 0.0519, 0.0458, 0.0406, 0.0358, 0.0316, 0.0277, 0.0243, 0.0213, 0.0185)

forecastedvalue = ft0 + pred5list

cor(forecastedvalue, real2018)^2

plot(modseason, which = 2)
shapiro.test(modseason$residuals)

qqnorm(mod200$residuals)
shapiro.test(mod200$residuals)

gt = predict(mod200, n.ahead = 1)
ft = fun1(t = 361, Month = 1)
ft+gt$pred[1]-(1.96*gt$se[1])
ft+gt$pred[1]+(1.96*gt$se[1])

co2.pred<-predict(mod200,n.ahead=12)
ft2 = fun1(t = c(361:372), Month = c(1:12))

plot(co2,xlim=c(0,375),ylim=c(350,420), xlab = "Months Since January 1988", ylab = "Levels of CO2 in the Atmosphere", main = "Level of CO2 in the Atmosphere Over Time", col = "blue")
lines((ft2+co2.pred$pred),col="red")
# Confidence interval
lines((ft2+co2.pred$pred)+1.96*co2.pred$se,col="orange",lty=3)
lines((ft2+co2.pred$pred)-1.96*co2.pred$se,col="orange",lty=3)

co2.pred<-predict(mod200,n.ahead=12)
ft3 = fun1(t = c(361:372), Month = c(1:12))

bestpred = ft3+co2.pred$pred

bestpred2 = c(407, bestpred)

sets = c(0,co2.pred$se)

plot(co2,xlim=c(325,375),ylim=c(380,420), xlab = "Months Since January 2015", ylab = "Levels of CO2 in the Atmosphere", main = "Level of CO2 in the Atmosphere Over Time", col = "blue")
#co2.pred<-predict(mod200,n.ahead=12)
lines(ts(bestpred2, start = 360),col="red")
# Confidence interval
lines(ts(bestpred2, start = 360)+1.96*sets,col="orange",lty=3)
lines(ts(bestpred2, start = 360)-1.96*sets,col="orange",lty=3)

layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot.ts(co2);acf(co2,main='Auto-Corr');pacf(co2,main='Part-Corr')

mod010 <- Arima(co2,order=c(0,1,0))
plot(mod010$residuals)
adf.test(mod010$residuals)

acf(mod010$residuals)
pacf(mod010$residuals)

mod210 <- Arima(co2,order=c(2,1,0))
acf(mod210$residuals)
pacf(mod210$residuals)

mod211 <- Arima(co2,order=c(2,1,1))
acf(mod211$residuals)
pacf(mod211$residuals)
coeftest(mod211)
summary(mod211)

mod311 <- Arima(co2,order=c(3,1,1))
acf(mod311$residuals)
pacf(mod311$residuals)
coeftest(mod311)
summary(mod311)

mod212 <- Arima(co2,order=c(2,1,2))
acf(mod212$residuals)
pacf(mod212$residuals)
coeftest(mod212)
summary(mod212)

auto.arima(co2)
bestmod <- Arima(co2,order=c(2,1,1))
acf(bestmod$residuals)

mod211010 <- Arima(co2,order=c(2,1,1),seasonal = list(order=c(1,1,2),period=12))
acf(mod211010$residuals)
pacf(mod211010$residuals)

season3fit <- ts(data.frame(co2, mod211010$fitted))
plot(season3fit,plot.type="s",lwd=3,col=c("red","purple"),ylab="CO2 levels",main="CO2 levels with SARIMA Model")

#Using auto.arima to check our work
co222 = ts(CO2Hawaii$CO2, frequency = 12)
asarima = auto.arima(co222, seasonal = T)

summary(asarima)

season4fit <- ts(data.frame(co2, asarima$fitted))
plot(season4fit,plot.type="s",lwd=3,col=c("red","purple"),ylab="CO2 levels",main="CO2 levels with Auto SARIMA Model")

pred6 <- forecast(mod211010, h=12)
cor(real2018, pred6$mean)^2

pred7 <- predict(mod211010, n.ahead=12)

bestpred = c(407, pred7$pred)

sets = c(0,pred7$se)

plot(co2,xlim=c(0,375),ylim=c(350,420), xlab = "Months Since January 1988", ylab = "Levels of CO2 in the Atmosphere", main = "Level of CO2 in the Atmosphere Over Time", col = "blue")
lines(ts(bestpred, start = 360),col="red")
# Confidence interval
lines(ts(bestpred, start = 360)+1.96*sets,col="orange",lty=3)
lines(ts(bestpred, start = 360)-1.96*sets,col="orange",lty=3)

qqnorm(mod211010$residuals)

pred6 <- forecast(mod211010, h=12)
cor(real2018, pred6$mean)^2

co2.pred<-predict(mod200,n.ahead=12)

ft3 = fun1(t = c(361:372), Month = c(1:12))

bestpred = ft3+co2.pred$pred

cor(real2018, bestpred)^2

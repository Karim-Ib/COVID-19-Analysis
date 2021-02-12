data_I = read.csv("C:\\Users\\Nutzer\\Documents\\untitled1\\time_series_covid19_confirmed_global.csv")
data_R = read.csv("C:\\Users\\Nutzer\\Documents\\untitled1\\time_series_covid19_recovered_global.csv")
head(data_R)

test.data_I = data_I[data_I$Country.Region == "Afghanistan", ]
test.data_R = data_R[data_R$Country.Region == "Afghanistan", ]



test.data_I = test.data_I[5:length(test.data_I)]
test.data_R = test.data_R[5:length(test.data_R)]
ind = which(test.data_I >= 100)[1]

test.data_I = as.numeric(test.data_I[ind:length(test.data_I)])
test.data_R = as.numeric(test.data_R[ind:length(test.data_R)])


write.csv(test.data_I, "C:\\Users\\Nutzer\\Documents\\untitled1\\afgahnistan.data.I.csv", row.names=FALSE, col.names=FALSE)
write.csv(test.data_R, "C:\\Users\\Nutzer\\Documents\\untitled1\\afgahnistan.data.R.csv", row.names=FALSE, col.names=FALSE)

plot(1:length(test.data_I), test.data_I)

#### Testdaten Italien


test.data_I = data_I[data_I$Country.Region == "Italy", ]
test.data_R = data_R[data_R$Country.Region == "Italy", ]

test.data_I = test.data_I[5:length(test.data_I)]
test.data_R = test.data_R[5:length(test.data_R)]
#plot(1:length(test.data_I),test.data_I)
ind_min = which(test.data_I >= 20)[1]
ind_max =65

png("C:\\Users\\Nutzer\\Documents\\untitled1\\italy_cases_data.png")
#plot(1:length(test.data_I),test.data_I, main="Covid cases in Sweden", xlab="Days",pch=16, ylab="Cases", cex=0.1, col="blue", type="p", yaxt="n")
plot(1:length(test.data_I),test.data_I, main="Covid cases in Italy", xlab="Days",pch=16, ylab="Cases", cex=0.1, col="blue", type="p")
#axis(2, at = c(1E5, 2E5, 3E5, 4E5, 5E5, 6E5), labels = c("100000", "200000", "300000", "400000", "500000", "600000"))
rect(xleft=c(ind_min, ind_min), xright=c(ind_max, ind_max), ytop=unlist(c(test.data_I[ind_max], ybottom=test.data_I[ind_max] )) ,c(0, 0),  lwd=2, border="red")
grid()
text(x=ind_min+10, y=test.data_I[ind_max]+50000, labels="Data before first lockdown", cex=0.7)
dev.off()

#test.data_I = as.numeric(test.data_I[ind:length(test.data_I)])
#test.data_R = as.numeric(test.data_R[ind:length(test.data_R)])
test.data_I = as.numeric(test.data_I[ind_min:ind_max])
test.data_R = as.numeric(test.data_R[ind_min:ind_max])


write.csv(test.data_I, "C:\\Users\\Nutzer\\Documents\\untitled1\\italy.data.I.csv", row.names=FALSE, col.names=FALSE)
write.csv(test.data_R, "C:\\Users\\Nutzer\\Documents\\untitled1\\italy.data.R.csv", row.names=FALSE, col.names=FALSE)

plot(1:length(test.data_I), test.data_I)


#### Testdaten Schweden


test.data_I = data_I[data_I$Country.Region == "Sweden", ]
test.data_R = data_R[data_R$Country.Region == "Sweden", ]


test.data_I = test.data_I[5:length(test.data_I)]
test.data_R = test.data_R[5:length(test.data_R)]
plot(1:length(test.data_I),test.data_I)
ind_min = which(test.data_I >= 25)[1]
ind_max = 100

png("C:\\Users\\Nutzer\\Documents\\untitled1\\sweden_cases_data.png")
plot(1:length(test.data_I),test.data_I, main="Covid cases in Sweden", xlab="Days",pch=16, ylab="Cases", cex=0.1, col="blue", type="p", yaxt="n")
axis(2, at = c(1E5, 2E5, 3E5, 4E5, 5E5, 6E5), labels = c("100000", "200000", "300000", "400000", "500000", "600000"))
rect(xleft=c(ind_min, ind_min), xright=c(ind_max, ind_max), ytop=unlist(c(test.data_I[ind_max], ybottom=test.data_I[ind_max] )) ,c(0, 0),  lwd=2, border="red")
grid()
text(x=ind_min, y=test.data_I[ind_max]+50000, labels="Data before first lockdown", cex=0.7)
dev.off()


#test.data_I = as.numeric(test.data_I[ind:length(test.data_I)])
#test.data_R = as.numeric(test.data_R[ind:length(test.data_R)])
test.data_I = as.numeric(test.data_I[ind_min:ind_max])
test.data_R = as.numeric(test.data_R[ind_min:ind_max])


write.csv(test.data_I, "C:\\Users\\Nutzer\\Documents\\untitled1\\sweden.data.I.csv", row.names=FALSE, col.names=FALSE)
write.csv(test.data_R, "C:\\Users\\Nutzer\\Documents\\untitled1\\sweden.data.R.csv", row.names=FALSE, col.names=FALSE)


#### Testdaten Spain


test.data_I = data_I[data_I$Country.Region == "Spain", ]
test.data_R = data_R[data_R$Country.Region == "Spain", ]


test.data_I = test.data_I[5:length(test.data_I)]
test.data_R = test.data_R[5:length(test.data_R)]

ind_min = which(test.data_I >= 20)[1]
ind_max = 80

plot(1:length(test.data_I),test.data_I)


png("C:\\Users\\Nutzer\\Documents\\untitled1\\spain_cases_data.png")
plot(1:length(test.data_I),test.data_I, main="Covid cases in Spain", xlab="Days",pch=16, ylab="Cases", cex=0.1, col="blue", type="p")
#axis(2, at = c(1E5, 2E5, 3E5, 4E5, 5E5, 6E5), labels = c("100000", "200000", "300000", "400000", "500000", "600000"))
rect(xleft=c(ind_min, ind_min), xright=c(ind_max, ind_max), ytop=unlist(c(test.data_I[ind_max], ybottom=test.data_I[ind_max] )) ,c(0, 0),  lwd=2, border="red")
grid()
text(x=ind_min, y=test.data_I[ind_max]+50000, labels="Data before first lockdown", cex=0.7)
dev.off()


#test.data_I = as.numeric(test.data_I[ind:length(test.data_I)])
#test.data_R = as.numeric(test.data_R[ind:length(test.data_R)])
test.data_I = as.numeric(test.data_I[ind_min:ind_max])
test.data_R = as.numeric(test.data_R[ind_min:ind_max])


write.csv(test.data_I, "C:\\Users\\Nutzer\\Documents\\untitled1\\spain.data.I.csv", row.names=FALSE, col.names=FALSE)
write.csv(test.data_R, "C:\\Users\\Nutzer\\Documents\\untitled1\\spain.data.R.csv", row.names=FALSE, col.names=FALSE)





#### Testdaten Deutschland


test.data_I = data_I[data_I$Country.Region == "Germany", ]
test.data_R = data_R[data_R$Country.Region == "Germany", ]


test.data_I = test.data_I[5:length(test.data_I)]
test.data_R = test.data_R[5:length(test.data_R)]
ind_min = which(test.data_I >= 20)[1]
ind_max = 80 #65


png("C:\\Users\\Nutzer\\Documents\\untitled1\\germany_cases_data.png")
plot(1:length(test.data_I),test.data_I, main="Covid cases in Germany", xlab="Days",pch=16, ylab="Cases", cex=0.1, col="blue", type="p")
rect(xleft=c(ind_min, ind_min), xright=c(ind_max, ind_max), ytop=unlist(c(test.data_I[ind_max], ybottom=test.data_I[ind_max] )) ,c(0, 0),  lwd=2, border="red")
grid()
text(x=ind_min, y=test.data_I[ind_max]+50000, labels="Data before first lockdown", cex=0.7)
dev.off()

#test.data_I = as.numeric(test.data_I[ind:length(test.data_I)])
#test.data_R = as.numeric(test.data_R[ind:length(test.data_R)])
test.data_I = as.numeric(test.data_I[ind_min:ind_max])
test.data_R = as.numeric(test.data_R[ind_min:ind_max])


write.csv(test.data_I, "C:\\Users\\Nutzer\\Documents\\untitled1\\germany.data.I.csv", row.names=FALSE, col.names=FALSE)
write.csv(test.data_R, "C:\\Users\\Nutzer\\Documents\\untitled1\\germany.data.R.csv", row.names=FALSE, col.names=FALSE)



### Intensiev Deutschland

germany.its =  read.csv("C:\\Users\\Nutzer\\Documents\\untitled1\\time_series_cvoid19_intensiev_germany.csv")
colnames(germany.its) = c("Datum", "Count")
head(germany.its)
count = germany.its$Count

ind.start = 170
ind.end = 270
germany.its$Datum[ind.start]
germany.its$Datum[ind.end]
datum.start = as.Date(germany.its$Datum[ind.start], format="%Y-%m-%d")
datum.end = as.Date(germany.its$Datum[ind.end], format="%Y-%m-%d")
#x.axis = seq(datum.start, datum.end, 1)

count = count[ind.start:ind.end]
plot(seq(datum.start, datum.end, 1), count)

write.csv(count, "C:\\Users\\Nutzer\\Documents\\untitled1\\germany.its.csv", row.names=FALSE, col.names=FALSE)



### 6th degree ICU (intensive care unit) regression for germany @2nd wave
d = 6
beta = c(2.36278271e+02, -2.17254582e+01,  4.30062690e+00, -2.60781616e-01, 6.86957012e-03, -7.44928865e-05,  2.82494897e-07)

x.axis = seq(datum.start, datum.end+30, 1)#datum.start + 0:150

X = matrix(rep(1, times=(d+1)*length(x.axis)), ncol=d+1)


for (i in 2:(d+1)){
	for (j in 1:length(x.axis)){
		X[j, i] = j^(i-1)
	}

}
 
X

y = X %*% beta

sum.max.icu = 24071


#x.axis[which(abs(y-sum.max.icu) < 800)]
png("C:\\Users\\Nutzer\\Documents\\untitled1\\germany_icu_prediction.png")


plot(x.axis, y, type="l", main="Prediction of ICU beds taken in Germany", xlab="Date", ylab="Beds", col="blue")
points(seq(datum.start, datum.end, 1), count, pch=16, cex=0.5, col="red")
legend("topleft", legend=c("Prediction", "Data"), col=c("blue", "red"), pch=c(16,16))
segments(0, sum.max.icu, x.axis[length(x.axis)]+10, sum.max.icu, col="black")
segments(x.axis[which(abs(y-sum.max.icu) < 800)], -5, x.axis[which(abs(y-sum.max.icu) < 800)], sum.max.icu, col="black", lty="dashed")
text(x.axis[which(abs(y-sum.max.icu) < 800)], 30, as.character(x.axis[which(abs(y-sum.max.icu) < 800)]), cex=0.8)
text(x.axis[60], 26000, "Total ICU capacity", cex=0.8)
grid()


dev.off()


###Correlation Deutschland Infected~Temperatur im Sommer

temp1 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata1.txt", header = TRUE, sep = "", dec = ".")
temp1 = temp1[2:length(temp1$TX), c(2,6) ]
temp1[,2] = as.numeric(temp1[,2])
temp1 = temp1[1:499,]
temp2 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata2.txt", header = TRUE, sep = "", dec = ".", na.strings=" ")
temp2 = temp2[2:length(temp2$TX),c(2,6)  ]
temp2[,2] = as.numeric(temp2[,2])
temp3 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata3.txt", header = TRUE, sep = "", dec = ".")
temp3 = temp3[2:length(temp3$TX), c(2,6) ]
temp3[,2] = as.numeric(temp3[,2])
temp4 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata4.txt", header = TRUE, sep = "", dec = ".")
temp4 = temp4[2:length(temp4$TX),c(2,6)  ]
temp4[,2] = as.numeric(temp4[,2])
temp5 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata5.txt", header = TRUE, sep = "", dec = ".")
temp5 = temp5[2:length(temp5$TX),c(2,6)  ]
temp5[,2] = as.numeric(temp5[,2])
temp6 = read.table("C:\\Users\\Nutzer\\Documents\\untitled1\\Wetterdata6.txt", header = TRUE, sep = "", dec = ".")
temp6 = temp6[2:length(temp6$TX),c(2,6)  ]
temp6[,2] = as.numeric(temp6[,2])



T_data = data.frame(Datum = as.Date(temp1$JJJJMMDD, format="%Y%m%d"), Temp=(temp1[,2] + temp2[,2] + temp3[,2] + temp4[,2] + temp5[,2] + temp6[,2]) / 6)
T_data = T_data[order(T_data$Datum, decreasing=FALSE),]




#T_data = data.frame(Datum = as.Date(temp1$JJJJMMDD, format="%Y%m%d"), Temp = as.numeric(temp$TM))
#T_data = T_data[order(T_data$Datum, decreasing=FALSE),]


test.data_I = data_I[data_I$Country.Region == "Germany", ]
test.data_I = test.data_I[5:length(test.data_I)]


plot(1:length(test.data_I), test.data_I)


ind_a = 270
ind_b = 370
test.data_I = test.data_I[(ind_a-1):ind_b]
test.data_dif = diff(as.numeric(test.data_I))
date_a = as.Date(names(test.data_I)[2], format="X%m.%d") -365
date_b = as.Date(names(test.data_I)[length(test.data_I)], format="X%m.%d")


ind_T_a = which(T_data$Datum == date_a)
ind_T_b = which(T_data$Datum == date_b)
#ind_T_b = ind_T_a - abs(ind_a-ind_b)

T_data = T_data[(ind_T_a):ind_T_b, ]







T_data$Temp = (T_data$Temp - mean(T_data$Temp))/sqrt(sd(T_data$Temp))
test.data_dif = (test.data_dif- mean(test.data_dif)) / sqrt(sd(test.data_dif))
reg_coef = summary(lm(as.numeric(test.data_dif)~T_data$Temp))$coefficients[c(1,2), 1]

cor = unclass(cor.test(T_data$Temp, as.numeric(test.data_dif), method="pearson"))$estimate

plot(T_data$Temp,as.numeric(test.data_dif),  type="n", xlab="Normalized Temperature", ylab="Normalized Cases", main="Correlation between new Covid cases and temperature")
points(T_data$Temp,as.numeric(test.data_dif), pch=16, col="red")
grid()
abline(reg_coef[1], reg_coef[2], col="blue", lwd=2)
abline(h=0)
abline(v=0)
text(4, 250, paste("Perason correlation =", round(cor,5)))



plot(1:length(test.data_I), test.data_I)
plot(1:length(test.data_dif), T_data$Temp)
x = 1:101
reg_coef = summary(lm(x~T_data$Temp))$coefficients[c(1,2), 1]
abline(reg_coef[1], reg_coef[2])



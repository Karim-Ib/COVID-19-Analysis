### Estimation of SIR beta&gamma parameters
### Using https://arxiv.org/pdf/1605.01931.pdf par. 3.3

#install.packages("deSolve")
library(deSolve)


Infected <- unclass(read.csv("C:\\Users\\Nutzer\\Documents\\untitled1\\test.data.csv"))$x
day <- 1:length(Infected)
N <- 38000000
init <- c(S = N-1, I = 1, R = 0)


SIR <- function(time, state, parameters){
	par <- as.list(c(state, parameters))
	with(par, { dS <- -beta * S * I / N 
			dI <- beta * S * I /N - gamma * I
			dR <- gamma * I
			list(c(dS, dI, dR))})
}

RSS.SIR <- function(parameters){
	names(parameters) <- c("beta", "gamma")
	out <- ode(y = init, times = day, func = SIR, parms = parameters)
	fit <- out[,3]
	RSS <- sum((Infected - fit)^2)
	return(RSS)
}

optim(c(0.001,0.4), RSS.SIR, method = "BFGS")








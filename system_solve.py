import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

### load data from RKI excel file
#R_data = pd.ExcelFile(r'/home/karim/Downloads/Nowcasting_Zahlen.xlsx')
#R_values = pd.read_excel(R_data,"Nowcast_R")
#R_values = pd.DataFrame(R_values, columns = ["Punktschätzer der Reproduktionszahl R"])
#R_dates = pd.read_excel(R_data,"Nowcast_R")
#R_dates = pd.DataFrame(R_dates, columns = ["Datum des Erkrankungsbeginns"])

### convert data into suitable numpy format + remove nan's
#R_values = pd.DataFrame.to_numpy(R_values)
#R_dates = pd.DataFrame.to_numpy(R_dates)
#R0 = np.empty(len(R_values))
#Dates = np.empty(len(R_dates))
#for i in range(len(R_values)):
    #R0[i]=R_values[i]
    #Dates[i] = i
#temp = 0
#for i in range(len(R0)):
    #if np.isnan(R0[temp]):
        #R0=np.delete(R0,temp)
    #else:
        #temp = temp+1

### set parameters
### Daten aus Deutschland vom 03.02.2021
gamma = 1/14
mu = 11.5/1000 ## sterberate pro 1000 jahr 2018
nu = 9.5/1000  ### geburtenrate pro 1000 jahr 2018 source https://www.laenderdaten.info/Europa/Deutschland/bevoelkerungswachstum.php
#R_avg = sum(R0)/len(R0)
R_avg = 1.3
#beta = gamma * R_avg
N = 83020000
I_in = 191000
R_in = 2000000
D_in = 58992
E_in = 300000
S_in_SEIR = N - I_in - R_in - E_in
S_in = N - I_in - R_in - D_in
#deaths = 58992
death_days = 23
mortality = D_in / (R_in * death_days)
#print(mortality)
incubation = 1 / 5
recovery = 1 - mortality
days = 365
dt = 10000
t = np.linspace(0,days, dt)

def solve_SRI(t, S_in, R_in, I_in, N, mu,nu,R_avg,gamma, days):

    #def dS(S,t,nu,mu,beta,N, I):
       #return  nu * N - beta * S * I / N- mu  * S;
    #def dI(I,t, mu, beta, gamma, N, S):
        #return beta * S * I / N - gamma * I - mu * I;
    #def dR(R,t, mu, gamma, I):
        #return gamma * I - mu * R;

    def system_dynamic(t, y, N, nu, mu, R_avg, gamma):
        beta = R_avg * (mu**2 + gamma * mu ) / nu
        S, R, I = y
        dydt = [nu * N - beta * S * I / N- mu * S, gamma * I - mu * R, beta * S * I / N - gamma * I - mu * I]
        return dydt;
    def system_static(t, y, N, R_avg, gamma):
        beta = R_avg * gamma
        S, R, I = y
        dydt = [ - beta * S * I / N, gamma * I , beta * S * I / N - gamma * I ]
        return dydt;
    y0 = [S_in, R_in, I_in]

    #SRI = odeint(system, y0, t, args=(N, nu, mu, beta, gamma))

    if days > 365 :
        #SRI = solve_ivp(system_dynamic, [0,days], y0, t_eval=t, args=[N, nu, mu, R_avg, gamma])
        SRI = solve_ivp(system_static, [0, days], y0, t_eval=t, args=[N, R_avg, gamma])
    else:
        SRI = solve_ivp(system_static, [0, days], y0, t_eval=t, args=[N, R_avg, gamma])

    #plt.plot(SRI.t, SRI.y[0], c="r", label="Nicht erkrankte Bevölkerung S(t)")
    #plt.plot(SRI.t, SRI.y[1], c="g", label="Immunisierte/Tote R(t)")
    #plt.plot(SRI.t, SRI.y[2], c="b", label="Infizierte I(t)")
    #plt.plot(SRI.t, SRI.y[0]+SRI.y[1]+SRI.y[2], c="yellow", label="N(t)" )
    #plt.xlabel("t in Tagen")
    #plt.ylabel("Anzahl der Personen")
    #if days > 365 :
        #plt.title("Dynamic SRI Model used for R0 = " + str(R_avg))
    #else:
        #plt.title("Static SRI Model used for R0 = " + str(R_avg))
    #plt.legend()
    #plt.show()

    return SRI;

def solve_SRID(t, S_in, R_in, I_in, D_in, N, Mortality, R_avg,gamma, days):
    mu = Mortality

    #def dS(S,t,nu,mu,beta,N, I):
       #return  nu * N - beta * S * I / N- mu  * S;
    #def dI(I,t, mu, beta, gamma, N, S):
        #return beta * S * I / N - gamma * I - mu * I;
    #def dR(R,t, mu, gamma, I):
        #return gamma * I - mu * R;

    def system_SRID(t, y, N, mu, R_avg, gamma):
        beta = R_avg * gamma
        S, R, I, D = y
        dydt = [- beta * S * I / N , gamma * I, beta * S * I / N  - gamma * I - mu * I, mu * I]
        return dydt;

    y0 = [S_in, R_in, I_in, D_in]

    SRID = solve_ivp(system_SRID, [0, days], y0, t_eval=t, args=[N, mu, R_avg, gamma])

    plt.plot(SRID.t, SRID.y[0], c="r", label="Nicht erkrankte Bevölkerung S(t)")
    plt.plot(SRID.t, SRID.y[1], c="g", label="Immunisierte R(t)")
    plt.plot(SRID.t, SRID.y[2], c="b", label="Infizierte I(t)")
    plt.plot(SRID.t, SRID.y[3], c="orange", label="Verstorbene D(t)")
    #plt.plot(SRID.t, SRID.y[0]+SRI.y[1]+SRI.y[2], c="yellow", label="N(t)" )
    plt.xlabel("t in Tagen")
    plt.ylabel("Anzahl der Personen")
    plt.title("SRID Model used for R0 = " + str(R_avg))
    plt.legend()
    plt.show()

    return;

def solve_SEIR(t, S_in, R_in, I_in, E_in, N, Incubation, R_avg, gamma, days):
    mu = Incubation

    #def dS(S,t,nu,mu,beta,N, I):
       #return  nu * N - beta * S * I / N- mu  * S;
    #def dI(I,t, mu, beta, gamma, N, S):
        #return beta * S * I / N - gamma * I - mu * I;
    #def dR(R,t, mu, gamma, I):
        #return gamma * I - mu * R;

    def system_SEIR(t, y, N, mu, R_avg, gamma):
        beta = R_avg * gamma
        S, R, E, I = y
        dydt = [- beta * S * I / N, gamma * I, beta * S * I / N - mu * E, mu * E - gamma * I]
        return dydt;

    y0 = [S_in, R_in, I_in, E_in]

    SEIR = solve_ivp(system_SEIR, [0, days], y0, t_eval=t, args=[N, mu, R_avg, gamma])

    plt.plot(SEIR.t, SEIR.y[0], c="r", label="Nicht erkrankte Bevölkerung S(t)")
    plt.plot(SEIR.t, SEIR.y[1], c="g", label="Immunisierte R(t)")
    plt.plot(SEIR.t, SEIR.y[2], c="b", label="Infizierte I(t)")
    plt.plot(SEIR.t, SEIR.y[3], c="orange", label="Exponiert E(t)")
    #plt.plot(SRID.t, SRID.y[0]+SRI.y[1]+SRI.y[2], c="yellow", label="N(t)" )
    plt.xlabel("t in Tagen")
    plt.ylabel("Anzahl der Personen")
    plt.title("SEIR Model used for R0 = " + str(R_avg))
    plt.legend()
    plt.show()

    return SEIR;

def solve_SEIRD(t, S_in, R_in, I_in, E_in, D_in, N, Incubation, Mortality, R_avg, gamma, days):
    delta = Incubation
    mu = Mortality

    #def dS(S,t,nu,mu,beta,N, I):
       #return  nu * N - beta * S * I / N- mu  * S;
    #def dI(I,t, mu, beta, gamma, N, S):
        #return beta * S * I / N - gamma * I - mu * I;
    #def dR(R,t, mu, gamma, I):
        #return gamma * I - mu * R;

    def system_SEIRD(t, y, N, delta, mu,  R_avg, gamma):
        beta = R_avg * gamma
        S, R, E, I, D = y
        dydt = [- beta * S * I / N, gamma * I, beta * S * I / N - delta * E, delta * E - (gamma + mu) * I, mu * I ]
        return dydt;

    y0 = [S_in, R_in, I_in, E_in, D_in]

    SEIRD = solve_ivp(system_SEIRD, [0, days], y0, t_eval=t, args=[N, delta, mu, R_avg, gamma])

    plt.plot(SEIRD.t, SEIRD.y[0], c="r", label="Susceptible S(t)")
    plt.plot(SEIRD.t, SEIRD.y[1], c="g", label="Recovered R(t)")
    plt.plot(SEIRD.t, SEIRD.y[2], c="b", label="Infected I(t)")
    plt.plot(SEIRD.t, SEIRD.y[3], c="orange", label="Exposed E(t)")
    plt.plot(SEIRD.t, SEIRD.y[4], c="purple", label="Deaths D(t)")
    #plt.plot(SEIRD.t, SEIRD.y[1] * mu, c="black", label="Test")
    #plt.plot(SRID.t, SRID.y[0]+SRI.y[1]+SRI.y[2], c="yellow", label="N(t)" )
    plt.xlabel("Days")
    plt.ylabel("Population")
    plt.grid()
    plt.title("SEIRD Model used for R0 = " + str(R_avg))
    plt.legend()
    plt.show()

    return;


#solve_SRI(t, N-(27*20+15), 15, 27*20, N, 0, 0, 2.383, 0.0272, days)
#solve_SRID(t, S_in, R_in, I_in, D_in, N, mortality, R_avg, gamma, days)
#solve_SEIR(t, S_in_SEIR, R_in, I_in, E_in, N, incubation, R_avg, gamma, days)
#solve_SEIRD(t, S_in_SEIR, R_in, I_in, E_in, D_in, N, incubation, mortality, R_avg, gamma, days)



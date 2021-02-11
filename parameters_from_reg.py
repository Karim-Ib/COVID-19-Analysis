import numpy as np
import numpy.linalg  as alg
import matplotlib.pyplot as plt
import scipy.stats as sci
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import pandas as pd
from system_solve import solve_SRI


data_italy_I = pd.read_csv("italy.data.I.csv")
data_italy_R = pd.read_csv("italy.data.R.csv")
#data_spain_I = pd.read_csv("spain.data.I.csv")
#data_spain_R = pd.read_csv("spain.data.R.csv")
data_spain_I = pd.read_csv("spain.data.I.csv")
data_spain_R = pd.read_csv("spain.data.R.csv")
data_germany_I = pd.read_csv("germany.data.I.csv")
data_germany_R = pd.read_csv("germany.data.R.csv")

data_italy_I = np.reshape(data_italy_I.to_numpy(), data_italy_I.shape[0])
data_italy_R = np.reshape(data_italy_R.to_numpy(), data_italy_R.shape[0])
data_spain_I = np.reshape(data_spain_I.to_numpy(), data_spain_I.shape[0])
data_spain_R = np.reshape(data_spain_R.to_numpy(), data_spain_R.shape[0])
data_germany_I = np.reshape(data_germany_I.to_numpy(), data_germany_I.shape[0])
data_germany_R = np.reshape(data_germany_R.to_numpy(), data_germany_R.shape[0])




N_italy = 60360000
N_spain = 46940000
N_germany = 83020000

data_italy_S = N_italy - (data_italy_I + data_italy_R)
data_spain_S =  N_spain - (data_spain_I + data_spain_R)
data_germany_S = N_germany - (data_germany_I + data_germany_R)

d = 3

k_italy = data_italy_I.shape[0]
k_spain = data_spain_I.shape[0]
k_germany = data_germany_I.shape[0]

x_italy = np.arange(0, k_italy, 1)
x_spain = np.arange(0, k_spain, 1)
x_germany = np.arange(0, k_germany, 1)

R_italy_reg = np.zeros(d+1)
R_spain_reg = np.zeros(d+1)
R_germany_reg = np.zeros(d+1)

S_italy_reg = np.zeros(d+1)
S_spain_reg = np.zeros(d+1)
S_germany_reg_reg = np.zeros(d+1)

dR_italy_reg = np.zeros(d)
dR_spain_reg = np.zeros(d)
dR_germany_reg = np.zeros(d)

dS_italy_reg = np.zeros(d)
dS_spain_reg = np.zeros(d)
dS_germany_reg_reg = np.zeros(d)

dS_italy = np.zeros(k_italy)
dS_spain = np.zeros(k_spain)
dS_germany = np.zeros(k_germany)

dR_italy = np.zeros(k_italy)
dR_spain = np.zeros(k_spain)
dR_germany = np.zeros(k_germany)

gamma_italy = 0
gamma_spain = 0
gamma_germany = 0

beta_italy = 0
beta_spain = 0
beta_germany = 0


def poly_deriv(coeff):
    deriv = np.zeros(len(coeff))
    #for i in range(len(deriv)):
        #deriv[i] = i * coeff[i]

    deriv = coeff * np.arange(0, len(coeff), 1)

    return deriv[1:len(deriv)];


def poly_reg(data, d):
    #data = np.reshape(data.to_numpy(), data.shape[0])
    X = np.empty((len(data), d + 1))
    X[:, 0] = np.ones(len(data))

    for k in range(1, d+1):
        for j in range(0, len(data)):
            X[j, k] = j**k
            if X[j, k] < 0 :
                print(j, k, X[j, k])


    beta = np.dot(np.dot(alg.inv(np.dot(np.transpose(X), X)), np.transpose(X)), data)

    return(beta);

R_italy_reg = poly_reg(data_italy_R, d)
R_spain_reg = poly_reg(data_spain_R, d)
R_germany_reg = poly_reg(data_germany_R, d)

S_italy_reg = poly_reg(data_italy_S, d)
S_spain_reg = poly_reg(data_spain_S, d)
S_germany_reg = poly_reg(data_germany_S, d)

dR_italy_reg = poly_deriv(R_italy_reg)
dR_spain_reg = poly_deriv(R_spain_reg)
dR_germany_reg = poly_deriv(R_germany_reg)

dS_italy_reg = poly_deriv(S_italy_reg)
dS_spain_reg = poly_deriv(S_spain_reg)
dS_germany_reg = poly_deriv(S_germany_reg)


#print(dR_italy_reg, R_italy_reg)

for i in range(k_italy):
    for dim in range(d-1):
        dS_italy[i] = dS_italy_reg[dim] * (x_italy[i] ** dim) + dS_italy[i]
        dR_italy[i] = dR_italy_reg[dim] * (x_italy[i] ** dim) + dR_italy[i]
#gamma_italy = gamma_italy + (dR_italy[i] / data_italy_I[i])
#beta_italy = beta_italy + (-dS_italy[i] * N_italy / (data_italy_I[i] * data_italy_S[i]))

#gamma_italy = gamma_italy / k_italy
#beta_italy = beta_italy / k_italy

gamma_italy = np.mean(dR_italy) / np.mean(data_italy_I)
beta_italy = np.mean(-dS_italy) * N_italy / (np.mean(data_italy_I) * np.mean(data_italy_S))

for i in range(k_spain):
    for dim in range(d):
        dS_spain[i] = dS_spain_reg[dim] * (x_spain[i] ** dim) + dS_spain[i]
        dR_spain[i] = dR_spain_reg[dim] * (x_spain[i] ** dim) + dR_spain[i]

gamma_spain = np.mean(dR_spain) / np.mean(data_spain_I)
beta_spain = np.mean(-dS_spain * N_spain / np.mean(data_spain_I * data_spain_S))

for i in range(k_germany):
    for dim in range(d):
        dS_germany[i] = dS_germany_reg[dim] * (x_germany[i]) ** dim + dS_germany[i]
        dR_germany[i] = dR_germany_reg[dim] * (x_germany[i]) ** dim + dR_germany[i]

gamma_germany = np.mean(dR_germany) / np.mean(data_germany_I)
beta_germany = np.mean(-dS_germany) * N_germany / (np.mean(data_germany_I) * np.mean(data_germany_S))


print(gamma_italy, beta_italy, beta_italy/gamma_italy)
print(gamma_spain, beta_spain, beta_spain/gamma_spain)
print(gamma_germany, beta_germany, beta_germany/gamma_germany)

gamma_italy = 0.015033198661431324
beta_italy = 0.21639420824132663

fun_ger = lambda mult: np.sqrt(sum((np.cumsum(solve_SRI(np.linspace(0,k_germany, k_germany), data_germany_S[0], data_germany_R[0], data_germany_I[0] * mult, N_germany, 0, 0, beta_germany / gamma_germany, gamma_germany, k_germany).y[2])-data_germany_I)**2))
fun_ita = lambda mult: np.sqrt(sum((np.cumsum(solve_SRI(np.linspace(0,k_italy, k_italy), data_italy_S[0], data_italy_R[0], data_italy_I[0]*mult, N_italy, 0, 0, beta_italy/gamma_italy, gamma_italy, k_italy).y[2])-data_italy_I)**2))
fun_spa = lambda mult: np.sqrt(sum((np.cumsum(solve_SRI(np.linspace(0,k_spain, k_spain), data_spain_S[0], data_spain_R[0], data_spain_I[0]*mult, N_spain, 0, 0, beta_spain/gamma_spain, gamma_spain, k_spain).y[2])-data_spain_I)**2))


factor_ger = minimize_scalar(fun_ger).x
factor_ita = minimize_scalar(fun_ita).x
factor_spa = minimize_scalar(fun_spa).x


#days = k_italy
#dt = days
#t = np.linspace(0,days, dt)


sir_spain = solve_SRI(np.linspace(0,k_spain, k_spain), data_spain_S[0], data_spain_R[0], data_spain_I[0]*factor_spa, N_spain, 0, 0, beta_spain/gamma_spain, gamma_spain, k_spain)
sir_germany = solve_SRI(np.linspace(0,k_germany, k_germany), data_germany_S[0], data_germany_R[0], data_germany_I[0]*factor_ger, N_germany, 0, 0, beta_germany/gamma_germany, gamma_germany, k_germany)
sir_italy = solve_SRI(np.linspace(0,k_italy, k_italy), data_italy_S[0], data_italy_R[0], data_italy_I[0]*factor_ita, N_italy, 0, 0, beta_italy/gamma_italy, gamma_italy, k_italy)

#np.savetxt("spain.SIR.I.txt", np.cumsum(test.y[2]))

plt.scatter(x_germany, data_germany_I, c="blue", s=2, label="Data")
plt.plot(sir_germany.t, np.cumsum(sir_germany.y[2]), c="red", label="SRI" )
plt.title("SRI compared to real data for the first wave in Germany")
plt.xlabel("Days")
plt.ylabel("Total cases")
plt.grid()
plt.legend()
plt.show()

plt.scatter(x_italy, data_italy_I, c="blue", s=2, label="Data")
plt.plot(sir_italy.t, np.cumsum(sir_italy.y[2]), c="red", label="SRI" )
plt.title("SRI compared to real data for the first wave in Italy")
plt.xlabel("Days")
plt.ylabel("Total cases")
plt.grid()
plt.legend()
plt.show()

plt.scatter(x_spain, data_spain_I, c="blue", s=2, label="Data")
plt.plot(sir_spain.t, np.cumsum(sir_spain.y[2]), c="red", label="SRI" )
plt.title("SRI compared to real data for the first wave in Spain")
plt.xlabel("Days")
plt.ylabel("Total cases")
plt.grid()
plt.legend()
plt.show()



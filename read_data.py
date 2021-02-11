import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from system_solve import solve_SEIR

#data_I = pd.read_csv("afgahnistan.data.I.csv")
#data_R = pd.read_csv("afgahnistan.data.R.csv")

data_I = pd.read_csv("italy.data.I.csv")
data_R = pd.read_csv("italy.data.R.csv")

data_I = np.reshape(data_I.to_numpy(), data_I.shape[0])
data_R = np.reshape(data_R.to_numpy(), data_R.shape[0])

k = len(data_I)
x = np.arange(0, k, 1)

N = 860360000
S = N - (data_I + data_R)

dS = np.empty(k - 1)
dR = np.empty(k - 1)

for i in range(k-1):
    dS[i] = S[i + 1] - S[i]
    dR[i] = data_R[i + 1] - data_R[i]

gamma = np.mean(dR / data_I[1:k])
beta = np.mean(-dS * N / (data_I[1:k] * S[1:k]))

print(gamma, beta, beta/gamma)


#dt = 10000
#t = np.linspace(0, 2*k, dt)
#E = 3 * data_I
#incubation = 1 / 8
#solution = solve_SRI(t, S[0], data_R[0], data_I[0], N, 0, 0, beta/gamma, gamma, 4*k)
#solve_SRI(t, S_in, R_in, I_in, N, mu, nu, R_avg, gamma, days)
#solution = solve_SEIR(t, S[0], data_R[0], 100*data_I[0], E[0], N, incubation, beta/gamma, gamma, 2*k)

#plt.scatter(x, data_I, s=2, c="r")
#plt.plot(solution.t, solution.y[2])
#plt.show()
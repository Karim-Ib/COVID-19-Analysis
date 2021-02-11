import numpy as np
import matplotlib.pyplot as plt
#import scipy.stats as sci
import scipy.optimize as sci
from scipy.integrate import solve_ivp
import pandas as pd

#from system_solve import solve_SIR, solve_SRID, solve_SEIR, solve_SEIRD


data = pd.read_csv("test.data.csv")
data = np.reshape(data.to_numpy(), data.shape[0])
x = np.arange(0,len(data), 1 )

cutoff_a = 10
cutoff_b = 100         ### arbitrary - want to catch the exponential relation at the beginning

lin_reg = np.polyfit(np.log(x[cutoff_a:cutoff_b]), np.log(data[cutoff_a:cutoff_b]), 1)


#plt.plot(np.log(x[cutoff_a:cutoff_b]), np.log(data[cutoff_a:cutoff_b]))
#plt.plot(np.log(x[cutoff_a:cutoff_b]), lin_reg[0]*np.log(x[cutoff_a:cutoff_b]) + lin_reg[1], c="r")
#plt.plot(x[cutoff_a:cutoff_b], (data[cutoff_a:cutoff_b]))
#plt.xlabel("time t")
#plt.ylabel("log(I)")
#plt.show()

#m = lin_reg[0]  ### m = beta - gamma
#print(m)
############################


def loss(beta, gamma, data, initial):
    size = len(data)

    def SIR(y):
        S = y[0]
        I = y[1]
        R = y[2]
        return [-beta*S*I, beta*S*I-gamma*I, gamma*I]
    print("test")
    solution = solve_ivp(SIR, [0, size], initial, t_eval=np.arange(0, size, 1), vectorized=True)
    return np.sqrt(np.mean((solution.y[1] - data) ** 2))
gamma = 1 / 30
N = 38000000
initial = [N-1, 1, 0]

optimal = sci.minimize_scalar(loss, args=(gamma, data, initial), method='Brent')


#optimal = sci.minimize(loss, [0.01, 0.01],args=(data, initial),method='L-BFGS-B',bounds=[(0.0001, 0.5), (0.0001, 0.5)])

beta, gamma = optimal.x




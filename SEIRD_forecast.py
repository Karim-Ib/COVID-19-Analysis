import scipy.stats as sci
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from system_solve import solve_SEIRD


gamma = 0.0364
R_1 = 0.8
R_2 = 2.3
R_3 = 5.4

N = 83020000
I_in = N * 0.03
R_in = N * 0.01
D_in = N * 0.0001
E_in = N * 0.05
S_in_SEIR = N - I_in - R_in - E_in
death_days = 23
mortality = D_in / (R_in * death_days)
incubation = 1 / 5
recovery = 1 - mortality
days = 365
dt = 10000
t = np.linspace(0,days, dt)

f1 = solve_SEIRD(t, S_in_SEIR, R_in, I_in, E_in, D_in, N, incubation, mortality, R_1, gamma, days)
f2 = solve_SEIRD(t, S_in_SEIR, R_in, I_in, E_in, D_in, N, incubation, mortality, R_2, gamma, days)
f3 = solve_SEIRD(t, S_in_SEIR, R_in, I_in, E_in, D_in, N, incubation, mortality, R_3, gamma, days)



import numpy as np
import numpy.linalg  as alg
import pandas as pd
import matplotlib.pyplot as plt

data_italy = pd.read_csv("italy.data.I.csv")
data_sweden = pd.read_csv("spain.data.I.csv")
data_germany = pd.read_csv("germany.data.I.csv")
data_germany_its = pd.read_csv("germany.its.csv")
d = 3

def poly_reg(data, d, name):
    data = np.reshape(data.to_numpy(), data.shape[0])
    x = np.arange(0, len(data), 1)
    X = np.empty((len(data), d + 1))
    X[:, 0] = np.ones(len(data))

    for k in range(1, d+1):
        for j in range(0, len(data)):
            X[j, k] = j**k
            if X[j, k] < 0 :
                print(j, k, X[j, k])


    beta_1 = np.dot(np.dot(alg.inv(np.dot(np.transpose(X), X)), np.transpose(X)), data)
    y_hat = np.dot(X, beta_1)

    #print(X[20, :])
    print(beta_1)

    MSE = ((data - y_hat)**2).mean()
    #print(MSE)
    RMSE = np.sqrt(MSE)
    #print(RMSE)

    R_2 = 1 - sum((data-y_hat)**2) / sum((data-np.mean(data))**2)
    #print(R_2)

    R_2_adj = round(R_2 - (1 - R_2) * d / (len(data) - d - 1), 4)
    #print(R_2_adj)

    temp = (data - np.dot(X, beta_1))
    sig = np.sqrt(np.dot(np.transpose(temp), temp) / (len(data) - d - 1))
    #print(sig)

    text_xpos = 30
    text_ypos = data[text_xpos-30]


    plt.plot(x, y_hat, c="g", label="Regression")
    plt.plot(x, y_hat - sig, c="r", linestyle="--", lw=0.7, label="Error")
    plt.plot(x, y_hat + sig, c="r", linestyle="--", lw=0.7)
    plt.scatter(x, data, c="b", s = 2, label="Recorded Data")
    plt.text(text_xpos, text_ypos, "R*^2 =" + str(R_2_adj))
    plt.grid()
    plt.title(str(d)+"th degree Polynomial regression for Covid cases in "+name )
    #plt.title(str(d) + "th degree Polynomial regression for ICU cases in " + name)
    plt.xlabel("Days")
    plt.ylabel("Cases")
    plt.legend()
    plt.show()

    return();

poly_reg(data_sweden, d, "Spain")


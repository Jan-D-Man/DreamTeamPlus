import numpy as np
import matplotlib.pyplot as plt
import math

K=0.56 
c_v=3686
ro=1081
V=40
sigma=0.472
T_0=1

t_adim=2*T_0*c_v*ro/(sigma*V**2)
l_adim=K*t_adim/(c_v*ro)

N_v = 100        # mida del vector (la matriu �s 100x100)
Tc = 36.5   
dX=0.02/l_adim
dt = 0.5*dX**2
gamma = dt/(dX)**2

Xmax = 100
tmax = 200

T01 = Tc
T00 = Tc
T0m1 = Tc
Tj1 = [Tc]*(int(tmax/dt)+2)
Ti = [Tc]*(int(Xmax/dX)+3)
Tdef = np.zeros((int(tmax/dt)+2,int(Xmax/dX)+3))
Tdef[:, 0] = Tj1
Tdef[0, :] = Ti


def tempEulExp(T0m1, T00, T01, dt, dX):
    t = 0
    j = 0
    while (t <= tmax):
        x = 0
        i = 1
        Tdef[:, int(Xmax/dX)+2] = Tc
        while (x < Xmax):
                T01= Tdef[j,i+1]
                T00 = Tdef[j, i]
                T0m1 = Tdef[j, i-1]
                T10 = dt*(((T01 - 2*T00 + T0m1)/(dX)**2) + 1) + T00
                x = x + dX
                i = i + 1
                Tdef[j+1, i-1] = T10
                
        j = j + 1
        t = t + dt
        Tdef[:, int(Xmax/dX)+2] = Tc
        Tdef[:, 0] = Tc
    return Tdef

def jacknicholson(T0m1, T00, T01,gamma):
    diagonal   = (1 + 2*gamma) * np.ones(N_v) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1 
    adalt  = (-gamma) * np.ones(N_v - 1) #�diagonal' de dalt i abaix 
    abaix  = (-gamma) * np.ones(N_v - 1)

    A = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #constru�m A

    diagonal   = (1 - 2*gamma) * np.ones(N_v) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1 
    adalt  = (gamma) * np.ones(N_v - 1) #�diagonal' de dalt i abaix 
    abaix  = (gamma) * np.ones(N_v - 1)

    B = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #constru�m B


    pass



print(jacknicholson(T0m1, T00, T01, gamma))
import numpy as np

K = 1
Tc = 36.5
dt = 1

Xmax = 10
tmax = 20
dX = 1
T01 = Tc
T00 = Tc
T0m1 = Tc
Tj1 = [Tc]*(int(tmax/dt)+2)
Ti = [Tc]*(int(Xmax/dX)+1)
Tdef = np.zeros((int(Xmax/dX)+1,int(tmax/dt)+2))
Tdef[0, :] = Tj1
Tdef[:, 0] = Ti


def tempEulExp(T0m1, T00, T01, dt, dX):
    t = 0
    j = 1
    while (t <= tmax):
        x = 0
        i = 0
        Tdef[:, int(tmax/dt)+1] = Tc
        while (x < Xmax):
                T01= Tdef[i+1,j]
                T00 = Tdef[i, j]
                T0m1 = Tdef[i-1, j]
                T10 = dt*(((T01 - 2*T00 + T0m1)/(dX)**2) + 1) + T00
                x = x + dX
                i = i + 1
                Tdef[i, j] = T10
                Tdef[0,j] = Tc
                
        j = j + 1
        t = t + dt
    Tdef[:, int(tmax/dt)+1] = Tc
    return Tdef

print(tempEulExp(T0m1, T00, T01, dt, dX))
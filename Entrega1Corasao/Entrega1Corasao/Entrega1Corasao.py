import numpy as np

K = 1
Tc = 36.5
dt = 0.1

Xmax = 10
tmax = 20
dX = 0.1
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

print(tempEulExp(T0m1, T00, T01, dt, dX))
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
DeltaX=0.02/l_adim
DeltaT = 0.5*DeltaX**2
alpha = DeltaT/(DeltaX)**2


for i in range(1,101):#va de 1 fins a 100

    b = np.zeros(N_v) #construeixo el vector b

    # j = 1 i j = N-2 en Python (posicions especials)
    b[0]    = Tc + DeltaT + alpha * Tc      # j = 1
    b[N_v-1]  = Tc + DeltaT + alpha * Tc      # j = N-1 Aqu� s� que puc dir que acaba en N_v-1 perqu� he dit que N_v �s 100

    # j = 2 fins a N-2  -> Python: 1 fins a N-2
    b[1:N_v-1] = Tc + DeltaT

    # Matriu tridiagonal A, com �s sim�tric la constru�m aix�
    diagonal   = (1 + 2*alpha) * np.ones(N_v) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1 
    adalt  = (-alpha) * np.ones(N_v - 1) #�diagonal' de dalt i abaix 
    abaix  = (-alpha) * np.ones(N_v - 1)

    A = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #constru�m A

    # Soluci� del sistema
    x = np.linalg.solve(A, b)

    inicial_final=np.zeros(N_v)
    inicial_final[0]=alpha*Tc
    inicial_final[-1]=alpha*Tc

    b=x+np.ones(N_v)*DeltaT+inicial_final

print(x)


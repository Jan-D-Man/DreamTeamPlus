# import numpy as np
# import matplotlib.pyplot as plt
# import math

# K=0.56
# c_v=3686
# ro=1081
# V=40
# sigma=0.472
# l_0=0.02
# tmax = 0.025
# t_0=(l_0**2)*c_v*ro/(K)
# T_0=t_0*sigma*V**2/(2*c_v*ro*0.02**2)

# N_v = 100        # mida del vector (la matriu �s 100x100)
# Tc = 36.5  
# Tc_norm=Tc/T_0

# dX=(0.02/100)/l_0
# dt = 0.5*dX**2
# gamma = dt/(dX)**2

# T01 = Tc_norm
# T00 = Tc_norm
# T0m1 = Tc_norm
# Tj1 = [Tc]*(int(tmax/dt)+2)
# Ti = [Tc]*(int(l_0/dX)+3)
# Tdef = np.zeros((int(tmax/dt)+2,int(l_0/dX)+3))
# Tdef[:, 0] = Tj1
# Tdef[0, :] = Ti


# def tempEulExp(T0m1, T00, T01, dt, dX):
#     t = 0
#     j = 0
#     while (t <= tmax):
#         x = 0
#         i = 1
#         Tdef[:, int(l_0/dX)+2] = Tc
#         while (x < l_0):
#                 T01= Tdef[j,i+1]
#                 T00 = Tdef[j, i]
#                 T0m1 = Tdef[j, i-1]
#                 T10 = dt*(((T01 - 2*T00 + T0m1)/(dX)**2) + 1) + T00
#                 x = x + dX
#                 i = i + 1
#                 Tdef[j+1, i-1] = T10
                
#         j = j + 1
#         t = t + dt
#         Tdef[:, int(l_0/dX)+2] = Tc
#         Tdef[:, 0] = Tc
#     return Tdef




# print(tempEulExp(T0m1, T00, T01, dt, dX))


import numpy as np
import matplotlib.pyplot as plt
import math

K=0.56
c_v=3686
ro=1081
V=40
sigma=0.472
l_0=0.02
t_0=(l_0**2)*c_v*ro/(K)
T_0=t_0*sigma*V**2/(2*c_v*ro*0.02**2)

N_v = 100        # mida del vector (la matriu �s 100x100)
Tc = 36.5  
Tc_norm=Tc/T_0

DeltaX=(0.02/100)/l_0
DeltaT = 0.49*DeltaX**2
alpha = DeltaT/(DeltaX)**2

#CAS 1

c = np.zeros(N_v-1) #construeixo el vector b
# j = 1 i j = N-2 en Python (posicions especials)
c[0]  = DeltaT + alpha * Tc_norm      # j = 1
c[-1]  =   DeltaT + alpha * Tc_norm     # j = N-1 Aqu� s� que puc dir que acaba en N_v-1 perqu� he dit que N_v �s 100

  # j = 2 fins a N-2  -> Python: 1 fins a N-2 no inclou el últim
c[1:-1] = DeltaT
    #Ax=b

diagonal   = (1 - 2*alpha) * np.ones(N_v-1) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1
adalt  = (alpha) * np.ones(N_v - 2) #�diagonal' de dalt i abaix
abaix = (alpha) * np.ones(N_v - 2)

B = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #construïm A

d=np.ones(N_v-1)*Tc_norm
temps = []
for i in range(0,int(0.025/DeltaT)):
    # Soluci� del sistema
    T=B@d + c
    temps.append(T)
    d=T

punts=np.linspace(0,0.02,99)

temps = np.array(temps)
plt.plot(punts,T*T_0)
plt.show()

plt.figure(figsize=(8,5))
plt.imshow(temps * T_0, 
           extent=[0, 0.02, 0, 500], 
           aspect='auto', 
           origin='lower', 
           cmap='hot')

plt.colorbar(label='Temperatura (°C)')
plt.xlabel('Posició (m)')
plt.ylabel('Passos de temps')
plt.title('Mapa de calor de l’evolució de temperatura')
plt.show()
#CAS 2

# DeltaT = 0.5*DeltaX**2
# alpha = DeltaT/(DeltaX)**2
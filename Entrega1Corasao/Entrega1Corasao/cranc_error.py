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
DeltaT = 0.5*DeltaX**2
alpha = DeltaT/(DeltaX)**2

#CAS 1

c = np.zeros(N_v-1) #construeixo el vector b
# j = 1 i j = N-2 en Python (posicions especials)
c[0]  = 2*DeltaT + 2*alpha * Tc_norm      # j = 1
c[-1]  =   2*DeltaT + 2*alpha * Tc_norm     # j = N-1 Aqu� s� que puc dir que acaba en N_v-1 perqu� he dit que N_v �s 100

  # j = 2 fins a N-2  -> Python: 1 fins a N-2 no inclou el últim
c[1:-1] = 2*DeltaT
    #Ax=b

# Matriu tridiagonal A, com �s sim�tric la constru�m aix�
diagonal_1   = 2*(1 + alpha) * np.ones(N_v-1) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1
adalt_1  = (-alpha) * np.ones(N_v - 2) #�diagonal' de dalt i abaix
abaix_1  = (-alpha) * np.ones(N_v - 2)

A = np.diag(diagonal_1) + np.diag(adalt_1, 1) + np.diag(abaix_1, -1) #construïm A

Ainv=np.linalg.inv(A)

diagonal_2   = 2*(1 - alpha) * np.ones(N_v-1) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1
adalt_2  = (alpha) * np.ones(N_v - 2) #�diagonal' de dalt i abaix
abaix_2 = (alpha) * np.ones(N_v - 2)

B = np.diag(diagonal_2) + np.diag(adalt_2, 1) + np.diag(abaix_2, -1) #construïm A

d=np.ones(N_v-1)*Tc_norm

temps = []

for i in range(0,500):
    # Soluci� del sistema
    T=Ainv@B@d+Ainv@c
    temps.append(T)
    d=T

print(Ainv@A)

punts=np.linspace(0,0.02,99)

plt.plot(punts,T*T_0)
plt.show()
temps = np.array(temps)

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

#SOLUCIÓ ANALÍTICA

t=0.025

def sol_anal(x_i):
    n = np.arange(1, 10**3)
    return Tc_norm + np.sum(
    (2/(n*np.pi)) * (1 - (-1)**n) *
    ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi**2 * n**2)) *
    np.sin(n*np.pi*x_i) ) 

T_anal = []

pos=np.linspace(0.01,0.99,99)

for posi in pos:
    T_anal.append(sol_anal(posi))

Error=[]

for i in range(99):
    Error.append(np.abs(T_anal[i]-T[i]))

print(Error)

plt.plot(pos,Error)
plt.show()
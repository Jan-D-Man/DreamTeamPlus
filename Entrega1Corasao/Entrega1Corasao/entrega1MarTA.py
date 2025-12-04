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
CASLIMIT = False
N_v = 100        # mida del vector (la matriu �s 100x100)
Tc = 36.5   
Tc_norm=Tc/T_0

DeltaX=(0.02/100)/l_0
DeltaT = DeltaX**2
alpha = DeltaT/(DeltaX)**2


#CAS 1 

b = np.zeros(N_v-1) #construeixo el vector b
# j = 1 i j = N-2 en Python (posicions especials)
b[0]    = Tc_norm + DeltaT + alpha * Tc_norm      # j = 1
b[-1]  = Tc_norm + DeltaT + alpha * Tc_norm      # j = N-1 Aqu� s� que puc dir que acaba en N_v-1 perqu� he dit que N_v �s 100

  # j = 2 fins a N-2  -> Python: 1 fins a N-2 no inclou el últim
b[1:-1] = Tc_norm + DeltaT
    #Ax=b

# Matriu tridiagonal A, com �s sim�tric la constru�m aix�
diagonal   = (1 + 2*alpha) * np.ones(N_v-1) #diagonal -> np.ones ens construeix un vector de dimensi� N_v amb tot 1 
adalt  = (-alpha) * np.ones(N_v - 2) #�diagonal' de dalt i abaix 
abaix  = (-alpha) * np.ones(N_v - 2)

A = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #construïm A

for i in range(0,250):#va de 1 fins a 100

    # Soluci� del sistema
    x = np.linalg.solve(A, b)

    inicial_final=np.zeros(N_v-1)
    inicial_final[0]=alpha*Tc_norm
    inicial_final[-1]=alpha*Tc_norm

    b=x+np.ones(N_v-1)*DeltaT+inicial_final
    if (CASLIMIT and x.max()*T_0 >= 50):
        break
    a = i


print("Valor de t màxim: ", a*DeltaT)
punts=np.linspace(0,0.02,99)
plt.plot(punts,x*T_0)
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
    Error.append(np.abs(T_anal[i]-x[i]))

print(Error)

plt.plot(pos,Error)
plt.show()


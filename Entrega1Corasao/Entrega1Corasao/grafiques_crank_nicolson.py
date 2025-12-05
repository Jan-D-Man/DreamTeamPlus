
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.animation import FuncAnimation
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter



K=0.56
c_v=3686
ro=1081
V=40
sigma=0.472
l_0=0.02
t_0=(l_0**2)*c_v*ro/(K)
T_0=t_0*sigma*V**2/(2*c_v*ro*0.02**2)
t=0.025

N_v = 99 # (la matriu es 99x99 i farem N_v de dimensió)
Tc = 36.5  
Tc_norm=Tc/T_0


def sol_anal(x_i): #definim la solucio analitica
    n = np.arange(1, 10**3)
    return Tc_norm + np.sum(
        (2/(n*np.pi)) * (1 - (-1)**n) *
        ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi**2 * n**2)) *np.sin(n*np.pi*x_i) )

T_anal = []
pos=np.linspace(0.01,0.99,99) #no hem agafat els termes de les C.C


for posi in pos: #Guardem totes les T analitiques per cada valor de x
    T_anal.append(sol_anal(posi))


Grafics=[]
Error_tres = [] 
pari = [0.5,1]

for par in pari:

    DeltaX=(0.02/100)/l_0
    DeltaT = par*DeltaX**2
    alpha = DeltaT/(DeltaX)**2

    
    c = np.zeros(N_v) #creem el vector independent del sistema matricial
    c[0]  = 2*DeltaT + 2*alpha * Tc_norm    #la primera i ultima tenen les C.C  
    c[-1]  =   2*DeltaT + 2*alpha * Tc_norm     
    c[1:-1] = 2*DeltaT
        
    #Creem la matriu que multiplica el vector de les T i+1
    diagonal_1   = 2*(1 + alpha) * np.ones(N_v) 
    adalt_1  = (-alpha) * np.ones(N_v - 1) 
    abaix_1  = (-alpha) * np.ones(N_v - 1)
    A = np.diag(diagonal_1) + np.diag(adalt_1, 1) + np.diag(abaix_1, -1) #construïm A

    Ainv=np.linalg.inv(A) #Calculem la inversa

    #Creem la matriu que multiplica el vector de les T i
    diagonal_2   = 2*(1 - alpha) * np.ones(N_v) 
    adalt_2  = (alpha) * np.ones(N_v - 1) 
    abaix_2 = (alpha) * np.ones(N_v - 1)
    B = np.diag(diagonal_2) + np.diag(adalt_2, 1) + np.diag(abaix_2, -1) 

    #El vector que conté T i (al inici valen Tc)
    d=np.ones(N_v)*Tc_norm

    temps = []

    for i in range(0,int(0.025/DeltaT)): #resolem el sistema
        T=Ainv@B@d+Ainv@c
        temps.append(T)
        d=T #Per cada iteració la T anterior canvia

    Grafics.append(T)

    Error=[]
    
    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-T[i]))

    Error_tres.append(Error)

fig, ax = plt.subplots(figsize=(8,5))

punts=np.linspace(0.0002,0.0198,99)
ax.plot(punts, Grafics[0] * T_0, label='$ΔT=0,5·ΔX^2$', linewidth=3, color='orange')
ax.plot(punts, Grafics[1] * T_0, label='$ΔT=ΔX^2$', linestyle='--',color='blue')
    

ax.set_xlabel('Posició (m)')
ax.set_ylabel('Temperatura (ºC)')
ax.set_xlim(0, 0.02)
ax.set_ylim(36, 54)
ax.grid(True,linestyle='--')

#fiquem els ticks
ax.tick_params(axis='both', which='both', direction='in', length=6)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

#Notació científia
class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax.xaxis.set_major_formatter(formatter)
ax.legend(loc='upper right')
plt.show()  

#Representem l'error
fig_err, ax_err = plt.subplots(figsize=(8,5))
ax_err.plot(pos, Error_tres[1], label='$ΔT=ΔX^2$', linewidth=3, color='orange')
ax_err.plot(pos, Error_tres[0], label='$ΔT=0,50·ΔX^2$', linestyle='--',color='blue')





ax_err.set_xlabel('Posició normalitzada')
ax_err.set_ylabel('Error')
ax_err.set_xlim(0, 1)
ax_err.set_ylim(0, 2.5*10**-6)
ax_err.grid(True,linestyle='--')

ax_err.tick_params(axis='both', which='both', direction='in', length=6)
ax_err.xaxis.set_ticks_position('both')
ax_err.yaxis.set_ticks_position('both')

class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%2.1f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax_err.yaxis.set_major_formatter(formatter)
ax_err.legend(loc='upper left')
plt.show()

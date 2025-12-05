
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

N_v = 99        # (la matriu es 99x99 i farem N_v de dimensió)
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

    b = np.zeros(N_v) #construeixo el vector que conté els valors de T anteriors
    b[0]    = Tc_norm + DeltaT + alpha * Tc_norm   
    b[-1]  = Tc_norm + DeltaT + alpha * Tc_norm   
    b[1:-1] = Tc_norm + DeltaT

    #Construeixo la matriu 99x99 per trobar T a un temps i
    diagonal   = (1 + 2*alpha) * np.ones(N_v) #La diagonal principal
    adalt  = (-alpha) * np.ones(N_v - 1) #la diagonal secundaria de dalt
    abaix  = (-alpha) * np.ones(N_v - 1) #la diagonal secundaria de abaix

    A = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #construïm la matriu

    inicial_final=np.zeros(N_v)
    inicial_final[0]=alpha*Tc_norm
    inicial_final[-1]=alpha*Tc_norm

    for i in range(0,int(0.025/DeltaT)): #resolc la equació fins a el temps que volem
        x = np.linalg.solve(A, b) #resol el sistema automàticament

        b=x+np.ones(N_v)*DeltaT+inicial_final #cada cop que iterem la temperatura anterior canvia
        
    Grafics.append(x)

    Error=[]
    
    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-x[i]))

    Error_tres.append(Error)

fig, ax = plt.subplots(figsize=(8,5))

punts=np.linspace(0.0002,0.0198,99)

#Fem el gràfic
ax.plot(punts, Grafics[0] * T_0, label='$ΔT=0,5·ΔX^2$', linewidth=3, color='orange')
ax.plot(punts, Grafics[1] * T_0, label='$ΔT=ΔX^2$', linestyle='--',color='blue')
    

ax.set_xlabel('Posició (m)')
ax.set_ylabel('Temperatura (ºC)')
ax.set_xlim(0, 0.02)
ax.set_ylim(36, 54)
ax.grid(True,linestyle='--')

# Fiquem ticks
ax.tick_params(axis='both', which='both', direction='in', length=6)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax.xaxis.set_major_formatter(formatter)
ax.legend(loc='upper right')
plt.show()  

#Sol analitica

#Escollim quin volem respresentar
fig_err, ax_err = plt.subplots(figsize=(8,5))
ax_err.plot(pos, Error_tres[1], label='$ΔT=ΔX^2$', linewidth=2)
ax_err.plot(pos, Error_tres[0], label='$ΔT=0,50·ΔX^2$', linewidth=2)



ax_err.set_xlabel('Posició normalitzada')
ax_err.set_ylabel('Error')
ax_err.set_xlim(0, 1)
ax_err.set_ylim(0, 0.16*10**-4)
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


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

N_v = 99        # (la matriu es 99x99 i farem N_v-1 de dimensió)
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



Error_tres = [] 
pari = [0.25, 0.49, 0.51]
Temperatures=[]

for par in pari:

    DeltaX=(0.02/100)/l_0
    DeltaT = par*DeltaX**2
    alpha = DeltaT/(DeltaX)**2

    c = np.zeros(N_v) #construeixo el vector del terme independent de la equació 
    c[1:-1] = DeltaT 
    c[0]  = DeltaT + alpha * Tc_norm  #El primer i últim vector tenen termes extra
    c[-1]  =   DeltaT + alpha * Tc_norm     
    
        
    #Construeixo la matriu 99x99 per trobar T a un temps i+1
    diagonal   = (1 - 2*alpha) * np.ones(N_v) #diagonal principal
    adalt  = (alpha) * np.ones(N_v - 1) #diagonal secundria de adalt
    abaix = (alpha) * np.ones(N_v - 1) #diagonal secundaria d'abaix
    B = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) 

    #Construeixo el vector que conté els valors de temperatura del temps anterior (per i=0 C.I.)
    d=np.ones(N_v)*Tc_norm

    temps = []
    for i in range(0,int(0.025/DeltaT)): #resolc la equació fins a el temps que volem
        T=B@d + c
        temps.append(T)
        d=T #per cada iteració fiquem la T del temps anterior

    Temperatures.append(T)

punts=np.linspace(0.0002,0.0198,99)   

fig, ax = plt.subplots(figsize=(8,5))

# Dibujamos la temperatura
ax.plot(punts, Temperatures[0]* T_0, label='$ΔT=0,25·ΔX^2$', linewidth=3, color='orange')
ax.plot(punts, Temperatures[1]* T_0, label='$ΔT=0,49·ΔX^2$', linestyle='--',color='blue')
# Configuramos etiquetas, límites y cuadrícula

ax.set_xlabel('Posició (m)')
ax.set_ylabel('Temperatura (ºC)')
ax.set_xlim(0, 0.02)
ax.set_ylim(36, 54)
ax.grid(True,linestyle='--')

# Configuramos las marcas de los ticks
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

#SOLUCIÓ ANALITICA
   
    Error=[]
    
    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-T[i]))

    Error_tres.append(Error)

fig_err, ax_err = plt.subplots(figsize=(8,5))
ax_err.plot(pos, Error_tres[2], label='$ΔT=0,25·ΔX^2$', linewidth=2)
ax_err.plot(pos, Error_tres[0], label='$ΔT=0,49·ΔX^2$', linewidth=2)
ax_err.plot(pos, Error_tres[1], label='$ΔT=0,51·ΔX^2$', linewidth=2)



# Configuramos etiquetas, límites y cuadrícula

ax_err.set_xlabel('Posició normalitzada')
ax_err.set_ylabel('Error')
ax_err.set_xlim(0, 1)
ax_err.set_ylim(0, 8*10**-6)
ax_err.grid(True,linestyle='--')

# Configuramos las marcas de los ticks
ax_err.tick_params(axis='both', which='both', direction='in', length=6)
ax_err.xaxis.set_ticks_position('both')
ax_err.yaxis.set_ticks_position('both')

class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax_err.yaxis.set_major_formatter(formatter)
# Mostramos la gráfica
ax_err.legend(loc='upper left')
plt.show()

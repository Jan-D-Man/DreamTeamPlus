from re import X
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
CAS_LIMIT = True #activar er el cas en que el model para a l'arribar a la temperatura límit (50ºC)


N_v = 99        # (la matriu es 99x99 i farem N_v-1 de dimensió)
Tc = 36.5  
Tc_norm=Tc/T_0


par=0.25
def euler_explicit(par):
   
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
        if CAS_LIMIT and np.any(T*T_0 >= 50):
            print(f"La simulació s'atura al pas de temps: {i*DeltaT:.4f} s perquè s'ha assolit la temperatura límit de 50ºC.")
            break
        d=T #per cada iteració fiquem la T del temps anterior
    
    fig, ax = plt.subplots(figsize=(8,5))
    punts=np.linspace(0.0002,0.0198,99)   
    plt.plot(punts,T*T_0) #Represento la temperatura sense normalitzar per cada
    plt.xlabel('Posició (m)')
    plt.ylabel('Temperatura (ºC)')

    class CustomScalarFormatter(ScalarFormatter):
        def _set_format(self):
            self.format = "%0.2f"

    formatter = CustomScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    formatter.set_scientific(True)

    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlim(0, 0.02)
    ax.set_ylim(36, 52)
    ax.grid(True,linestyle='--')

    # Configuramos las marcas de los ticks
    ax.tick_params(axis='both', which='both', direction='in', length=6)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.legend(loc='upper right')

    plt.show()

  
euler_explicit(par)


    
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.animation import FuncAnimation


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
        d=T #per cada iteració fiquem la T del temps anterior

    punts=np.linspace(0,0.02,99)   
    plt.plot(punts,T*T_0) #Represento la temperatura sense normalitzar per cada
    plt.xlabel('Posició (m)')
    plt.ylabel('Passos de temps')
    plt.show()

    temps = np.array(temps)
    plt.figure(figsize=(8,5)) #Represento el diagrama de calor
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

    fig, ax = plt.subplots()

    # línia inicial (temps 0)
    linea, = ax.plot(punts, temps[0] * T_0)

    ax.set_xlim(0, 0.02)
    ax.set_ylim((temps*T_0).min()*0.9, (temps*T_0).max()*1.1)
    ax.set_xlabel("Posició (m)")
    ax.set_ylabel("Temperatura (°C)")

    def init():
        """Configura el primer frame"""
        linea.set_ydata(temps[0] * T_0)
        ax.set_title("t = 0.0 s")
        return linea,

    def update(frame):
        """Actualitza la línia a cada pas de temps"""
        T_actual = temps[frame] * T_0
        linea.set_ydata(T_actual)
        t_real = frame * DeltaT * t_0   # si vols temps físic
        ax.set_title(f"t = {t_real:.2f} s")
        return linea,

    ani = FuncAnimation(
        fig,
        update,
        frames=len(temps),
        init_func=init,
        blit=True,
        interval=50   # ms entre frames (canvia-ho per fer-la més ràpida/lenta)
    )

    plt.show()
    

    def sol_anal(x_i): #definim la solucio analitica
        n = np.arange(1, 10**3)
        return Tc_norm + np.sum(
        (2/(n*np.pi)) * (1 - (-1)**n) *
        ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi**2 * n**2)) *np.sin(n*np.pi*x_i) )

    T_anal = []
    pos=np.linspace(0.01,0.99,99)

    for posi in pos: #Guardem totes les T analitiques per cada valor de x
        T_anal.append(sol_anal(posi))

    Error=[]
    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-T[i])) 

    print(Error)

    plt.plot(pos,Error)
    plt.ylabel(Error)
    plt.xlabel(pos)
    plt.show()

par_1=1

def euler_implicit(par_1):
    DeltaX=(0.02/100)/l_0
    DeltaT = par_1*DeltaX**2
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

        b=x+np.ones(N_v-1)*DeltaT+inicial_final #cada cop que iterem la temperatura anterior canvia

    punts=np.linspace(0,0.02,99)
    plt.plot(punts,x*T_0) 
    plt.show()

   

    #SOLUCIÓ ANALÍTICA

    t=0.025

    def sol_anal(x_i):
        n = np.arange(1, 10**3)
        return Tc_norm + np.sum(
        (2/(n*np.pi)) * (1 - (-1)**n) *
        ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi * n**2)) *
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


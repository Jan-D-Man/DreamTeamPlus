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

def sol_anal(x_i): #definim la solucio analitica
    n = np.arange(1, 10**3)
    return Tc_norm + np.sum(
        (2/(n*np.pi)) * (1 - (-1)**n) *
        ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi**2 * n**2)) *np.sin(n*np.pi*x_i) )

T_anal = []
pos=np.linspace(0.01,0.99,99) #no hem agafat els termes de les C.C

for posi in pos: #Guardem totes les T analitiques per cada valor de x
    T_anal.append(sol_anal(posi))



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

    punts=np.linspace(0.0002,0.0198,99)   
    plt.plot(punts,T*T_0) #Represento la temperatura sense normalitzar per cada
    plt.xlabel('Posició (m)')
    plt.ylabel('Temperatura (ºC)')
    plt.title('Euler explicit')
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

    fig2, ax2 = plt.subplots(figsize=(8, 3))

# Rang dinàmic més sensible
    all_vals = temps * T_0
    vmin = np.percentile(all_vals, 1)
    vmax = np.percentile(all_vals, 99)

    # Primera fila (temps = 0)
    data0 = all_vals[0][np.newaxis, :]

    # Heatmap inicial
    img = ax2.imshow(
        data0,
        extent=[0, 0.02, 0, 1],
        aspect='auto',
        origin='lower',
        cmap='turbo',
        vmin=vmin,
        vmax=vmax
    )

    # Colorbar
    cbar = plt.colorbar(img, ax=ax2, label='Temperatura (°C)')
    ax2.set_xlabel('Posició (m)')
    ax2.set_yticks([])

    # --- LÍNIES DE FRONTERA DEL CENTRE (0.5 cm = 0.005 m) ---
    x_center = 0.01
    half = 0.005 / 2
    x_left = x_center - half       # 0.0075
    x_right = x_center + half      # 0.0125

    # Línies discontínues
    ax2.axvline(x=x_left, color='white', linestyle='--', linewidth=1.4)
    ax2.axvline(x=x_right, color='white', linestyle='--', linewidth=1.4)

    # --- FRANJA SEMITRANSPARENT ENTRE LES DUES LÍNIES ---
    ax2.axvspan(x_left, x_right, color='white', alpha=0.15)  
    # alpha = 0.15 dóna un ressalt suau sense tapar el mapa de calor

    ax2.set_title('t = 0.000 s')

    # Funció d'actualització
    def update_heat(frame):
        data = all_vals[frame][np.newaxis, :]
        img.set_data(data)
        t_real = frame * DeltaT
        ax2.set_title(f't = {t_real:.3f} s')
        return img,
    frames_to_show = range(0, len(temps), 8)
    ani2 = FuncAnimation(
        fig2,
        update_heat,
        frames=frames_to_show,
        interval=1,   # més ràpid
        blit=False
    )

    plt.show()
    
    #SOLUCIÓ ANALITICA
    Error=[]
    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-T[i])) 

    print(Error)

    plt.plot(pos,Error)
    plt.ylabel('Error')
    plt.xlabel('Posició normalitzada')
    plt.title('Error euler explicit ')
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

        b=x+np.ones(N_v)*DeltaT+inicial_final #cada cop que iterem la temperatura anterior canvia

    punts=np.linspace(0.0002,0.0198,99)
    plt.plot(punts,x*T_0) 
    plt.xlabel('Posició (m)')
    plt.ylabel('Temperatura (ºC)')
    plt.show()

   

    #SOLUCIÓ ANALÍTICA
    Error=[]

    for i in range(99): #Calculem l'error
        Error.append(np.abs(T_anal[i]-x[i]))

    plt.plot(pos,Error)
    plt.xlabel('Posició normalitzada')
    plt.ylabel('Error')
    plt.title('Error euler implícit')
    plt.show()


per_2=0.5

def crank_nicolson(per_2):

    DeltaX=(0.02/100)/l_0
    DeltaT = per_2*DeltaX**2
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


    punts=np.linspace(0.0002,0.0198,99)
    plt.plot(punts,T*T_0)
    plt.xlabel('Posició (m)')
    plt.ylabel('Temperatura (ºC)')
    plt.title('Cranck-Nicolson')
    plt.show()

    temps = np.array(temps)

    plt.figure(figsize=(8,5))
    plt.imshow(temps * T_0, 
               extent=[0.0002, 0.0198, 0, int(0.025/DeltaT)], 
               aspect='auto', 
               origin='lower', 
               cmap='hot')

    plt.colorbar(label='Temperatura (°C)')
    plt.xlabel('Posició (m)')
    plt.ylabel('Passos de temps')
    plt.title('Mapa de calor de l’evolució de temperatura')
    plt.show()

    #SOLUCIÓ ANALÍTICA
    Error=[]

    for i in range(99):
        Error.append(np.abs(T_anal[i]-T[i]))

    plt.plot(pos,Error)
    plt.xlabel('Posició normalitzada')
    plt.ylabel('Error')
    plt.title('Error Cranck-Nicolson')
    plt.show()

euler_explicit(par)
euler_implicit(par_1)
crank_nicolson(per_2)

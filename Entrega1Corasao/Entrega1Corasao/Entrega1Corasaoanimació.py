
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

N_v = 100        # mida del vector (la matriu �s 100x100)
Tc = 36.5  
Tc_norm=Tc/T_0

DeltaX=(0.02/100)/l_0
DeltaT = 0.25*DeltaX**2
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
t=0.025

def sol_anal(x_i):
    n = np.arange(1, 10**3)
    return Tc_norm + np.sum(
    (2/(n*np.pi)) * (1 - (-1)**n) *
    ((1 - np.exp(-n**2 * np.pi**2 * t)) / (np.pi**2 * n**2)) *np.sin(n*np.pi*x_i) )

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


#CAS 2

# DeltaT = 0.5*DeltaX**2
# alpha = DeltaT/(DeltaX)**2
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
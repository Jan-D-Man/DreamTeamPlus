
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
t=0.025


#ANIMACIÓ

# DeltaT = 0.5*DeltaX**2
# alpha = DeltaT/(DeltaX)**2
fig2, ax2 = plt.subplots(figsize=(8, 3))

all_vals = temps * T_0
vmin = np.percentile(all_vals, 1)
vmax = np.percentile(all_vals, 99)

data0 = all_vals[0][np.newaxis, :]

img = ax2.imshow(
    data0,
    extent=[0, 0.02, 0, 1],
    aspect='auto',
    origin='lower',
    cmap='turbo',
    vmin=vmin,
    vmax=vmax
)

cbar = plt.colorbar(img, ax=ax2, label='Temperatura (°C)')
ax2.set_xlabel('Posició (m)')
ax2.set_yticks([])
x_center = 0.01
half = 0.005 / 2
x_left = x_center - half       # 0.0075
x_right = x_center + half      # 0.0125

# Línies discontínues
ax2.axvline(x=x_left, color='white', linestyle='--', linewidth=1.4)
ax2.axvline(x=x_right, color='white', linestyle='--', linewidth=1.4)

ax2.axvspan(x_left, x_right, color='white', alpha=0.15)  

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

class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax2.xaxis.set_major_formatter(formatter)

plt.show()
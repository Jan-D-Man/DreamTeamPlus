import numpy as np
import matplotlib.pyplot as plt
import math

N_v = 100        # mida del vector (exemple)
Tc = 36.5      
DeltaT = 0.00025    
alpha = 1/0.00025     

for i in range(1,101):#va de 1 fins a 100

    b = np.zeros(N_v) #construeixo el vector b

    # j = 1 i j = N-2 en Python (posicions especials)
    b[0]    = Tc + DeltaT + alpha * Tc      # j = 1
    b[N_v-1]  = Tc + DeltaT + alpha * Tc      # j = N-1 Aquí sí que puc dir que acaba en N_v-1 perquè he dit que N_v és 100

    # j = 2 fins a N-2  -> Python: 1 fins a N-2
    b[1:N_v-1] = Tc + DeltaT

    # Matriu tridiagonal A, com és simètric la construïm així
    diagonal   = (1 + 2*alpha) * np.ones(N_v) #diagonal -> np.ones ens construeix un vector de dimensió N_v amb tot 1 
    adalt  = (-alpha) * np.ones(N_v - 1) #´diagonal' de dalt i abaix 
    abaix  = (-alpha) * np.ones(N_v - 1)

    A = np.diag(diagonal) + np.diag(adalt, 1) + np.diag(abaix, -1) #construïm A

    # Solució del sistema
    x = np.linalg.solve(A, b)

    b=x+np.ones(N_v)*(DeltaT+i*DeltaT+alpha*Tc)

print(x)


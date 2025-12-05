# DreamTeamPlus
DreamTeam + Jan D Man

S’ha d’instal·lar la llibreria Numpy i Mathplotlib per a poder visualitzar el nostre codi.
Al Github tenim diversos documents:
- Enterga definitiva (document on està el codi principal)
- Gràfiques euler explícit
- Gràfiques euler implícit
- Gràfiques crank-nicolson
- Solució del problema amb euler explícit
- Animacions
  
Entrega definitiva: aquí està el codi principal. Està dividit en tres parts definides: euler explícit, euler implícit
i crank-nicolson. Per a fer-los aparèixer (tant l’evolució de la temperatura en funció de la posició com l’error) s’ha de
cridar-los al final del codi amb els seg¨uents noms:

- euler explicit(par)
- euler implicit(par 1)
- crank nicolson(per 2)
  
On fica par, par1 i par2 ficar el valors que multiplica (∆˜x)2 que volguem veure representats. A l’Euler explícit
també aparèix l’animació del diagrama de la transmisió del calor entre les plaques.

Gràfiques euler explícit, implícit i crank-nicolson: Aquests tres codis els hem utilitzat per a presentar amb
la forma correcta els gràfics que estan al document ’latex’.

En els tres hem utilitzat bucles per a representar els gràfics amb els diversos paràmetres (el valors que multipliquen
(∆˜x)2). El gràfic Temperatura en funció de Posició del euler explícit l’hem hagut de separar en dos ja que el valor de
0.51 divergia. Per a veure un gràfic o l’altre treure ficar # a les dues primeres línies o només a la última (gràfic que
divergeix):

- ax.plot(punts, Temperatures[0]* T0, label=’∆T = 0.25∆X2’, linewidth=3, color=’orange’)
- ax.plot(punts, Temperatures[1]* T0, label=’∆T = 0, 49∆X2’, linestyle=’–’,color=’blue’)
- ax.plot(punts, Temperatures[1]* T0, label=’∆T = 0, 51∆X2’, color=’blue’)
  
Els límits per al primer cas han de ser plt.ylim(36,53) i pel segon plt.ylim(-1000,1000). No cal tocar el límit en l’eix x.
Amb l’error també passa el mateix. Seleccionar amb # les dues primeres o la última:

- ax err.plot(pos, Error tres[2], label=’∆T = 0, 25∆X2’, linewidth=2)
- ax err.plot(pos, Error tres[0], label=’∆T = 0, 49∆X2’, linewidth=2)
- ax err.plot(pos, Error tres[1], label=’∆T = 0, 51∆X2’, linewidth=2)
  
Si seleccionem les primeres deixar el plt.ylim(0,8*10**-6) i si volem representar el 0.51 posar plt.ylim(0,1.6). No
canviar el límit en l’eix x.

Solució del problema amb euler explícit: per a veure el cas solucionat cal anar o bé a l'arxiu solució_problema o bé al de entrega definitva i canviar la variable CAS_LIMIT a True.

A l’apartat d’animacions estan:
- L'animació de la transmisió del calor entre les plaques per l’euler explícit


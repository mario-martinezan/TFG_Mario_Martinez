import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
import csv

universe = mda.Universe("NP_HSA-prot_autopsf.pdb", ["MD_s1.dcd"])

proteina = universe.select_atoms("protein")
radio = []
tiempo = []
frames = []
trayectoria = universe.trajectory
print("frames = ", len(trayectoria))

for i in universe.trajectory:
    frame = universe.trajectory.frame
    t = universe.trajectory.time
    r = universe.atoms.radius_of_gyration()
    tiempo.append(t)
    radio.append(r/10) #en nanometros
    frames.append(frame)

with open("datos_radius+NP.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Frame', 'Tiempo (ps)', 'Radio (nm)']) 
    for f, t, r in zip(frames, tiempo, radio):
        writer.writerow([f, t, r]) 

print("todo ok")

eqs= sum(radio[2500:])/len(radio[2500:])
eq= [eqs]*len(tiempo[2500:])
print("Radius at equilibrium (nm) = ", eqs)

plt.figure(figsize=(8, 6))
plt.plot(tiempo, radio, label="Radius of gyration (conformational change) (nm)", color="orange", linewidth=1)
plt.plot(tiempo[2500:], eq, label="Radius at equilibrium (nm)", color="black", linewidth=0.8)
plt.title("Radius of gyration (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.legend()
plt.savefig("f_radius+NP")

###############################################

r2= eqs
r1 = np.linspace(3,9.4)
numero_max = []

def n_max (x): #x es r2/r1
    alpha = math.asin(x/(1+x))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax


for i in r1:
    numero_max.append(n_max(r2/i))

x_sim = np.array([3, 6, 9]) 
y_sim = np.array([3, 10, 30]) 
coeficientes_sim = np.polyfit(x_sim, y_sim, 2)  # Ajuste de curva
polinomio_sim = np.poly1d(coeficientes_sim)  # Crear la función del polinomio
curva_sim = polinomio_sim(r1) 

x_ex = np.array([4, 5.7, 9.4]) 
y_ex = np.array([4, 12, 40]) 
coeficientes_ex = np.polyfit(x_ex, y_ex, 2)  # Ajuste de curva
polinomio_ex = np.poly1d(coeficientes_ex)  # Crear la función del polinomio
curva_ex = polinomio_ex(r1) 

plt.figure(figsize=(8, 6))
plt.title("Maximun number of proteins in the surface", fontsize=16, pad=15)
plt.xlabel("Nanoparticle radius (nm)", fontsize=14, labelpad=10)
plt.ylabel("Nº of proteins", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

plt.plot(r1, numero_max, color="orange", linewidth=1)
plt.scatter(3,n_max(r2/3), color= "orange")
plt.scatter(6,n_max(r2/6), color= "orange")
plt.scatter(9,n_max(r2/9), color= "orange", label="Model (with conformational change)")

plt.plot(r1, curva_sim, color="green", linewidth=1)
plt.scatter(3,3, color= "green")
plt.scatter(6,10, color= "green")
plt.scatter(9,30, color= "green", label="Simulation")

plt.plot(r1, curva_ex, color="blue", linewidth=1)
plt.scatter(4,4, color= "blue")
plt.scatter(5.7,12, color= "blue")
plt.scatter(9.4,40, color= "blue", label="Experimental")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("f_nmax_radius+NP")
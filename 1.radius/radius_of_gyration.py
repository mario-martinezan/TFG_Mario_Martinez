import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np


universe = mda.Universe("HSA_autopsf.pdb", ["simulacio_HSA_part1.dcd","simulacio_HSA_part2.dcd"])
#universe = mda.Universe("HSA_autopsf.pdb", "simulacio_HSA_part1.dcd")

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
    radio.append(r/10)
    frames.append(frame)

eqs= sum(radio[3500:])/len(radio[3500:])
eq= [eqs]*len(tiempo[3500:])
print("Radius at equilibrium (nm) = ", eqs)

plt.figure(figsize=(8, 6))
plt.plot(tiempo, radio, label="Radius of gyration (nm)", color="red", linewidth=1)
plt.plot(tiempo[3500:], eq, label="Radius at equilibrium (nm)", color="black", linewidth=0.8)
plt.title("Radius of gyration (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.legend()
plt.savefig("Radius_of_gyration")

###############################################

r2= 2.71
r1 = np.linspace(3,9.4)
numero_max = []

def n_max (x):
    alpha = math.asin(x/ (1+x))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax


for i in r1:
    numero_max.append(n_max(r2/i))

x_sim = np.array([3, 6, 9]) 
y_sim = np.array([3, 10, 30]) 
coeficientes = np.polyfit(x_sim, y_sim, 2)  
polinomio = np.poly1d(coeficientes) 
curva_sim = polinomio(r1) 

x_ex = np.array([4, 5.7, 9.4]) 
y_ex = np.array([4, 12, 40]) 
coeficientes_ex = np.polyfit(x_ex, y_ex, 2)  
polinomio_ex = np.poly1d(coeficientes_ex) 
curva_ex = polinomio_ex(r1) 

plt.figure(figsize=(8, 6))
plt.title("Maximun number of proteins in the surface", fontsize=16, pad=15)
plt.xlabel("Nanoparticle radius (nm)", fontsize=14, labelpad=10)
plt.ylabel("Nº of proteins", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

plt.plot(r1, numero_max, color="red", linewidth=1)
plt.scatter(3,n_max(r2/3), color= "red")
plt.scatter(6,n_max(r2/6), color= "red")
plt.scatter(9,n_max(r2/9), color= "red", label="Model")

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
plt.savefig("n_proteins")

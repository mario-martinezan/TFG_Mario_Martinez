import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
import csv

universe = mda.Universe("HSA_autopsf.pdb", ["simulacio_HSA_part1.dcd","simulacio_HSA_part2.dcd"])
#universe = mda.Universe("HSA_autopsf.pdb", "simulacio_HSA_part1.dcd")

protein = universe.select_atoms("protein")
trayectoria = universe.trajectory

distancia = []
tiempo = []
frames = []
print("frames = ", len(trayectoria))

"""
for i in universe.trajectory:
    frame = universe.trajectory.frame
    t = universe.trajectory.time
    positions = protein.positions
    distances = pdist(positions)
    max_distance = np.max(distances) #no es la radio, es la distancia
    tiempo.append(t)
    distancia.append(max_distance/10)
    frames.append(frame)


with open("datos_max_distance.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Frame', 'Tiempo (ps)', 'Distancia máxima (nm)']) 
    for f, t, d in zip(frames, tiempo, distancia):
        writer.writerow([f, t, d]) 

print("dooc guardado")

"""
with open('datos_max_distance.csv', mode='r') as file:
    reader = csv.reader(file)
    header = next(reader)  
    for row in reader:
        tiempo.append(float(row[1]))  
        distancia.append(float(row[2]))  
        frames.append(int(row[0]))  

print("todo ok")


eqs= sum(distancia[3500:])/len(distancia[3500:])
print("Max. distance at equilibrium (nm) = ", eqs)


plt.figure(figsize=(8, 6))
plt.plot(tiempo, distancia, label="Max. distance between two atoms (nm)", color="purple", linewidth=1)
plt.plot(tiempo[3500:], [eqs]*len(tiempo[3500:]), label="Max. distance at equilibrium (nm)", color="black", linewidth=0.8)
plt.title("Max. distance between two atoms in HSA (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Max. Distance (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
#plt.legend()
plt.savefig("f_max_distance")

#################################3

r2= eqs/2 #ahora si radio
r1 = np.linspace(3,9.4,500)
numero_max = []

#la x es el radio de la NP
def n_max (x):
    alpha = math.asin(r2/(x+r2)) #cambio respecto a q sea una esfera, ahora es "mas plana"
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax


for i in r1:
    numero_max.append(n_max(i))

x_sim = np.array([3, 6, 9]) 
y_sim = np.array([3, 10, 30]) 
coeficientes = np.polyfit(x_sim, y_sim, 2)  # Ajuste de curva
polinomio = np.poly1d(coeficientes)  # Crear la función del polinomio
curva_sim = polinomio(r1) 

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

plt.plot(r1, curva_sim, color="green", linewidth=1)
plt.scatter(3,3, color= "green")
plt.scatter(6,10, color= "green")
plt.scatter(9,30, color= "green", label="Simulation")

plt.plot(r1, curva_ex, color="blue", linewidth=1)
plt.scatter(4,4, color= "blue")
plt.scatter(5.7,12, color= "blue")
plt.scatter(9.4,40, color= "blue", label="Experimental")

plt.plot(r1, numero_max, color="purple", linewidth=1)
plt.scatter(3,n_max(3), color= "purple")
plt.scatter(9,n_max(9), color= "purple", label="Model (max. distance)")
plt.scatter(6, n_max(6), color= "purple")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("f_nmax_max_distance")

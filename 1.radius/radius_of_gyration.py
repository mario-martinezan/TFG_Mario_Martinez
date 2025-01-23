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

plt.figure(figsize=(10, 6))
plt.plot(tiempo, radio, label="Radius of gyration (nm)", color="red", linewidth=1)
plt.title("Radius of gyration (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("Radius_of_gyration")

###############################################

r2= 2.71
r1 = np.linspace(1,10)
numero_max = []

def n_max (x):
    alpha = math.asin(x/ (1+x))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax


for i in r1:
    numero_max.append(n_max(r2/i))

plt.figure(figsize=(10, 6))
plt.plot(r1, numero_max, color="red", linewidth=1)
plt.title("Maximun number of proteins in the surface", fontsize=16, pad=15)
plt.xlabel("Nanoparticle radius (nm)", fontsize=14, labelpad=10)
plt.ylabel("Nº of proteins", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.scatter(3,3, color= "green")
plt.scatter(6,10, color= "green")
plt.scatter(9,30, color= "green", label="Simulation")
plt.scatter(3,n_max(r2/3), color= "red")
plt.scatter(6,n_max(r2/6), color= "red")
plt.scatter(9,n_max(r2/9), color= "red", label="Model")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("n_proteins")

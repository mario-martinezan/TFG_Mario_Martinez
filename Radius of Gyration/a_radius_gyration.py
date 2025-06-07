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
plt.plot(tiempo, radio, label="Radius of gyration (nm)", color="orange", linewidth=1)
plt.plot(tiempo[3500:], eq, label="Radius at equilibrium (nm)", color="black", linewidth=1.8)
#plt.title("Radius of gyration (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.legend()
plt.savefig("f_radius_of_gyration")
plt.savefig("f_radius_of_gyration.eps", format="eps")

###############################################

r2= eqs
r1 = np.linspace(3/2,9.4/2)
numero_max = []


def n_max (x):
    alpha = math.asin(x/ (1+x))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax

con_masa = []
def n_max_mas (m, r1):
    #print(m)
    #m1 = (m / (1.66054 * 10**(-21)) ) #* (6.02214076*10**23)
    #print(m1)
    m1 = m
    r2 = 0.066 * (m1**(1/3))
    #print(r2)
    alpha = math.asin(r2/(r2+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    #print(nmax)
    return nmax

masa = proteina.total_mass() # en umas?
for i in r1:
    numero_max.append(n_max(r2/i))
    con_masa.append(n_max_mas(masa, i))

x_sim = np.array([3/2, 6/2, 9/2]) 
y_sim = np.array([3, 10, 30]) 
coeficientes = np.polyfit(x_sim, y_sim, 2)  # Ajuste de curva
polinomio = np.poly1d(coeficientes)  # Crear la función del polinomio
curva_sim = polinomio(r1) 

x_ex = np.array([4/2, 5.7/2, 9.4/2]) 
y_ex = np.array([4, 12, 40]) 
coeficientes_ex = np.polyfit(x_ex, y_ex, 2)  # Ajuste de curva
polinomio_ex = np.poly1d(coeficientes_ex)  # Crear la función del polinomio
curva_ex = polinomio_ex(r1) 

plt.figure(figsize=(8, 6))
#plt.title("Maximun number of proteins in the surface", fontsize=16, pad=15)
plt.xlabel("Nanoparticle radius (nm)", fontsize=14, labelpad=10)
plt.ylabel("# of proteins", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

plt.plot(r1, numero_max, color="orange", linewidth=2)
plt.scatter(3/2,n_max(r2/(3/2)), color= "orange", linewidth=4)
plt.scatter(6/2,n_max(r2/(6/2)), color= "orange", linewidth=4)
plt.scatter(9/2,n_max(r2/(9/2)), color= "orange", linewidth=4, label="Radius of gyration Model")

plt.plot(r1, con_masa, color="red", linewidth=2)
plt.scatter(3/2,n_max_mas(masa, (3/2)), color= "red", linewidth=4)
plt.scatter(6/2,n_max_mas(masa, (6/2)), color= "red", linewidth=4)
plt.scatter(9/2,n_max_mas(masa, (9/2)), color= "red", linewidth=4, label="Geometric Model")

plt.plot(r1, curva_sim, color="purple", linewidth=2)
plt.scatter(3/2,3, color= "purple", linewidth=4)
plt.scatter(6/2,10, color= "purple", linewidth=4)
plt.scatter(9/2,30, color= "purple", linewidth=4, label="MD Simulation")

plt.plot(r1, curva_ex, color="black", linewidth=2)
plt.scatter(4/2,4, color= "black", linewidth=4)
plt.scatter(5.7/2,12, color= "black", linewidth=4)
plt.scatter(9.4/2,40, color= "black", linewidth=4, label="Experimental")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig("f_radius_gy")
plt.savefig("f_radius_gy.eps", format="eps")


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
masa = protein.total_mass() # en umas?
print("masa ", masa)

r_gran = []
r_mid = []
r_peq = []
tiempo = []
frames = []
print("frames = ", len(trayectoria))


for i in universe.trajectory:
    frame = universe.trajectory.frame
    t = universe.trajectory.time
    tiempo.append(t)
    frames.append(frame)
    I_tensor = protein.moment_of_inertia()
    I_values = np.linalg.eigvals(I_tensor) /100 # Tres valores correspondientes a los ejes principales EN NANOMETROS
    radio1 = np.sqrt(I_values[0]*5/(2*masa) ) # Calcular los radios de inercia como si fuese una espera solida a partir de los ejer principales
    radio2 = np.sqrt(I_values[1]*5/(2*masa) ) 
    radio3 = np.sqrt(I_values[2]*5/(2*masa) ) 
    if radio1>radio2 and radio1>radio3:
        r_gran.append(radio1)
        if radio2>radio3:
            r_mid.append(radio2)
            r_peq.append(radio3)
        else:
            r_mid.append(radio3)
            r_peq.append(radio2)
    elif radio2>radio1 and radio2>radio3:
        r_gran.append(radio2)
        if radio1>radio3:
            r_mid.append(radio1)
            r_peq.append(radio3)
        else:
            r_mid.append(radio3)
            r_peq.append(radio1)
    else:
        r_gran.append(radio3)
        if radio1>radio2:
            r_mid.append(radio1)
            r_peq.append(radio2)
        else:
            r_mid.append(radio2)
            r_peq.append(radio1)
    


with open("datos_inertia_radius.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Frame', 'Tiempo (ps)', 'Radio 1 (nm)', 'Radio 2 (nm)', 'Radio 3 (nm)']) 
    for f, t, r1, r2, r3 in zip(frames, tiempo, r_gran, r_mid, r_peq):
        writer.writerow([f, t, r1, r2, r3]) 

print("todo ok")


"""
with open('datos_inertia_radius.csv', mode='r') as file:
    reader = csv.reader(file)
    header = next(reader)  
    for row in reader:
        frames.append(int(row[0])) 
        tiempo.append(float(row[1]))  
        r_gran.append(float(row[2])) 
        r_mid.append(float(row[3])) 
        r_peq.append(float(row[4]))  
 

print("todo ok")
"""

eqs1= sum(r_gran[3500:])/len(r_gran[3500:])
eqs2= sum(r_mid[3500:])/len(r_mid[3500:])
eqs3= sum(r_peq[3500:])/len(r_peq[3500:])
print("Inertia radius at equilibrium (nm) = [ ", eqs1, " , ", eqs2, " , ", eqs3, " ]")


plt.figure(figsize=(8, 6))
plt.plot(tiempo, r_gran, label=" Radius of inertia 1 (nm)", color="cyan", linewidth=1)
plt.plot(tiempo, r_mid, label=" Radius of inertia 2 (nm)", color="magenta", linewidth=1)
plt.plot(tiempo, r_peq, label=" Radius of inertia 3 (nm)", color="yellow", linewidth=1)

plt.plot(tiempo[3500:], [eqs1]*len(tiempo[3500:]), color="black", linewidth=0.8)
plt.plot(tiempo[3500:], [eqs2]*len(tiempo[3500:]), color="black", linewidth=0.8)
plt.plot(tiempo[3500:], [eqs3]*len(tiempo[3500:]), color="black", linewidth=0.8)

plt.title("Radius of inertia  (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
#plt.legend()
plt.savefig("f_inertia_radius")


#################################3


r1 = np.linspace(3,9.4)
nradio1 = []
nradio2 = []
nradio3 = []

def n_max (r2,r1):
    alpha = math.asin(r2/(r2+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax


for i in np.linspace(3,9.4):
    nradio1.append(n_max(eqs1,i))
    nradio2.append(n_max(eqs2,i))
    nradio3.append(n_max(eqs3,i))

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

plt.plot(r1, nradio1, color="cyan", linewidth=1)
plt.scatter(3,n_max(eqs1,3), color= "cyan")
plt.scatter(6,n_max(eqs1,6), color= "cyan")
plt.scatter(9,n_max(eqs1,9), color= "cyan", label="Model (Radius 1 of inercia)")

plt.plot(r1, nradio2, color="magenta", linewidth=1)
plt.scatter(3,n_max(eqs2,3), color= "magenta")
plt.scatter(6,n_max(eqs2,6), color= "magenta")
plt.scatter(9,n_max(eqs2,9), color= "magenta", label="Model (Radius 2 of inercia)")

plt.plot(r1, nradio3, color="yellow", linewidth=1)
plt.scatter(3,n_max(eqs3,3), color= "yellow")
plt.scatter(6,n_max(eqs3,6), color= "yellow")
plt.scatter(9,n_max(eqs3,9), color= "yellow", label="Model (Radius 3 of inercia)")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("f_nmax_inertia_radius")
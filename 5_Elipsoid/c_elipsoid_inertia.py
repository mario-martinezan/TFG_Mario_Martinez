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


#calculos momentos inercia elipsoide
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
    if I_values[0]>I_values[1] and I_values[0]>I_values[2]:
        Ivalue0 = I_values[0]
        if I_values[1]>I_values[2]:
            Ivalue1 = I_values[1]
            Ivalue2 = I_values[2]
        else:
            Ivalue1 = I_values[2]
            Ivalue2 = I_values[1]
    elif I_values[1]>I_values[0] and I_values[1]>I_values[2]:
        Ivalue0 = I_values[1]
        if I_values[0]>I_values[2]:
            Ivalue1 = I_values[0]
            Ivalue2 = I_values[2]
        else:
            Ivalue1 = I_values[2]
            Ivalue2 = I_values[0]
    else:
        Ivalue0 = I_values[2]
        if I_values[0]>I_values[1]:
            Ivalue1 = I_values[0]
            Ivalue2 = I_values[1]
        else:
            Ivalue1 = I_values[1]
            Ivalue2 = I_values[0]
    radio1 = np.sqrt(Ivalue0*5/(2*masa) ) 
    radio2 = np.sqrt((Ivalue1*5/masa)-radio1**2) 
    radio3 = np.sqrt((Ivalue2*5/masa)-radio1**2) 
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
    


with open("datos_inertia_radius_elip.csv", mode='w', newline='') as file:
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

eqs1= sum(r_gran[3500:])/len(r_gran[3500:]) #radio mayor
eqs2= sum(r_mid[3500:])/len(r_mid[3500:])  #radio medio
eqs3= sum(r_peq[3500:])/len(r_peq[3500:])
print("Inertia radius at equilibrium (nm) = [ ", eqs1, " , ", eqs2, " , ", eqs3, " ]")

plt.figure(figsize=(8, 6))
plt.plot(tiempo, r_gran, label=" Radius of inertia 1 (nm)", color="cyan", linewidth=1)
plt.plot(tiempo, r_mid, label=" Radius of inertia 2 (nm)", color="magenta", linewidth=1)
plt.plot(tiempo, r_peq, label=" Radius of inertia 3 (nm)", color="yellow", linewidth=1)

plt.plot(tiempo[3500:], [eqs1]*len(tiempo[3500:]), color="black", linewidth=0.8)
plt.plot(tiempo[3500:], [eqs2]*len(tiempo[3500:]), color="black", linewidth=0.8)
plt.plot(tiempo[3500:], [eqs3]*len(tiempo[3500:]), color="black", linewidth=0.8)

plt.title("Radius of inertia for an elipsoid (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
#plt.legend()
plt.savefig("f_inertia_radius_elip")


#################################3
r_granc = []
r_midc = []
r_peqc = []
tiempoc = []
framesc = []
print("frames = ", len(trayectoria))


with open('datos_inertia_radius.csv', mode='r') as file:
    reader = csv.reader(file)
    header = next(reader)  
    for row in reader:
        framesc.append(int(row[0])) 
        tiempoc.append(float(row[1]))  
        r_granc.append(float(row[2])) 
        r_midc.append(float(row[3])) 
        r_peqc.append(float(row[4])) 

eqs1c= sum(r_granc[3500:])/len(r_granc[3500:]) #radio mayor
eqs2c= sum(r_midc[3500:])/len(r_midc[3500:])  #radio medio
eqs3c= sum(r_peqc[3500:])/len(r_peqc[3500:])


#################33
r1 = np.linspace(3,9.4)

def n_max_elip (r1, a,b):
    alpha = math.asin((b)/(a+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax

nradio = []
nradioc = []

for i in r1:
    nradio.append(n_max_elip(i, eqs2,eqs1))
    nradioc.append(n_max_elip(i, eqs2c,eqs1c))

#def inversa(r):
 #   return 1 - (r - 3) / (9.4 - 3)



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

plt.plot(r1, nradio, color="chocolate", linewidth=1)
plt.scatter(3,n_max_elip(3, eqs2, eqs1), color= "chocolate")
plt.scatter(6,n_max_elip(6, eqs2, eqs1), color= "chocolate")
plt.scatter(9,n_max_elip(9, eqs2, eqs1), color= "chocolate", label="Model (Elipsoid with I elip.)")

plt.plot(r1, nradioc, color="pink", linewidth=1)
plt.scatter(3,n_max_elip(3, eqs2c, eqs1c), color= "pink")
plt.scatter(6,n_max_elip(6, eqs2c, eqs1c), color= "pink")
plt.scatter(9,n_max_elip(9, eqs2c, eqs1c), color= "pink", label="Model (Elipsoid with I esfe.)")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("f_nmax_elipsoid")
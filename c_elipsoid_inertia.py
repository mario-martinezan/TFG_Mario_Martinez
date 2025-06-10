import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
import csv
#from matplotlib.lines import Line2D

universe = mda.Universe("HSA_autopsf.pdb", ["simulacio_HSA_part1.dcd","simulacio_HSA_part2.dcd"])
#universe = mda.Universe("HSA_autopsf.pdb", "simulacio_HSA_part1.dcd")

protein = universe.select_atoms("protein")
trayectoria = universe.trajectory
masa = protein.total_mass() # en umas?
print("masa ", masa)


#calculos momentos inercia elipsoide
r_1 = []
r_2 = []
r_3 = []
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
    I_values.sort()  #ordena de menor a mayor
    radio1 = np.sqrt(I_values[0]*5/(2*masa) )           # el radio peuqeño, repetido b=c
    radio2 = np.sqrt((I_values[1]*5/masa)-radio1**2)    # el radio grande
    radio0 = np.sqrt((I_values[2]*5/masa)-radio1**2)    # el radio grande pero otra opción
    radios = []
    radios.append(radio1)
    radios.append(radio2)
    radios.append(radio0)
    radios.sort()
    r_1.append(radios[1])  # esto es b, radio "grande" del elipsodie
    r_2.append(radios[0])   # esto es a, radio "pequeño" del elipsodie
    r_3.append(radios[2])
    #r_1.append(radio2)  # esto es b, radio "grande" del elipsodie
    #r_2.append(radio1)   # esto es a, radio "pequeño" del elipsodie
    #r_3.append(radio0)

    


with open("datos_inertia_radius_elip.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Frame', 'Tiempo (ps)', 'Radio 1 (nm)', 'Radio 2 (nm)', 'Radio 3 (nm)']) 
    for f, t, r1, r2, r3 in zip(frames, tiempo, r_1, r_2, r_3):
        writer.writerow([f, t, r1, r2, r3]) 

print("todo ok")


"""
with open('datos_inertia_radius.csv', mode='r') as file:
    reader = csv.reader(file)
    header = next(reader)  
    for row in reader:
        frames.append(int(row[0])) 
        tiempo.append(float(row[1]))  
        r_1.append(float(row[2])) 
        r_2.append(float(row[3])) 
        r_3.append(float(row[4]))  
 

print("todo ok")
"""

eqs1= sum(r_1[3500:])/len(r_1[3500:]) 
eqs2= sum(r_2[3500:])/len(r_2[3500:]) 
eqs3= sum(r_3[3500:])/len(r_3[3500:])
print("Inertia radius at equilibrium (nm) = [ ", eqs1, " , ", eqs2, " , ", eqs3, " ]")

plt.figure(figsize=(8, 6))
plt.plot(tiempo, r_3, label=" Radius of inertia 1 (nm)", color="springgreen", linewidth=1)
plt.plot(tiempo, r_1, label=" Radius of inertia 2 (nm)", color="green", linewidth=1)
plt.plot(tiempo, r_2, label=" Radius of inertia 3 (nm)", color="seagreen", linewidth=1)


plt.plot(tiempo[3500:], [eqs1]*len(tiempo[3500:]), color="black", linewidth=1, label="Radius at equilibrium (nm)")
plt.plot(tiempo[3500:], [eqs2]*len(tiempo[3500:]), color="black", linewidth=1)
plt.plot(tiempo[3500:], [eqs3]*len(tiempo[3500:]), color="black", linewidth=1)

#plt.title("Radius of inertia for an elipsoid (nm)", fontsize=16, pad=15)
plt.xlabel("Time (ps)", fontsize=14, labelpad=10)
plt.ylabel("Radius (nm)", fontsize=14, labelpad=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.legend(fontsize=13, loc='right', bbox_to_anchor=(1, 0.61))
plt.savefig("f_inertia_radius_elip")




#################################3
r_1c = []
r_2c = []
r_3c = []
tiempoc = []
framesc = []
print("frames = ", len(trayectoria))


with open('datos_inertia_radius.csv', mode='r') as file:
    reader = csv.reader(file)
    header = next(reader)  
    for row in reader:
        framesc.append(int(row[0])) 
        tiempoc.append(float(row[1]))  
        r_1c.append(float(row[2])) 
        r_2c.append(float(row[3])) 
        r_3c.append(float(row[4])) 

eqs1c= sum(r_1c[3500:])/len(r_1c[3500:]) #radio mayor
eqs2c= sum(r_2c[3500:])/len(r_2c[3500:])  #radio medio
eqs3c= sum(r_3c[3500:])/len(r_3c[3500:])


#################33
s3= 3/2
s6= 6/2
s9= 9/2
r1 = np.linspace(s3,4.7,50)

def n_max_elip (r1, a,b):
    alpha = math.asin((b)/(a+r1))
    nmax1 = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    nmax = nmax1
    return nmax

nradio = []
nradioc = []

for i in r1:
    nradio.append(n_max_elip(i, eqs2,eqs1))
    nradioc.append(n_max_elip(i, eqs2c,eqs1c))

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

masa = protein.total_mass() # en umas?
for i in r1:
    con_masa.append(n_max_mas(masa, i))

x_sim = np.array([s3, s6, s9]) 
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

plt.plot(r1, curva_sim, color="purple", linewidth=2)
plt.scatter(s3,3, color= "purple", linewidth=4)
plt.scatter(s6,10, color= "purple", linewidth=4)
plt.scatter(s9,30, color= "purple", linewidth=4, label="MD Simulation")

plt.plot(r1, curva_ex, color="black", linewidth=2)
plt.scatter(4/2,4, color= "black", linewidth=4)
plt.scatter(5.7/2,12, color= "black", linewidth=4)
plt.scatter(9.4/2,40, color= "black", linewidth=4, label="Experimental")

plt.plot(r1, nradio, color="green", linewidth=2)
plt.scatter(s3,n_max_elip(s3, eqs2, eqs1), color= "green", linewidth=4)
plt.scatter(s6,n_max_elip(s6, eqs2, eqs1), color= "green", linewidth=4)
plt.scatter(s9,n_max_elip(s9, eqs2, eqs1), color= "green", linewidth=4, label="Elipsoidal Model")

plt.plot(r1, con_masa, color="red", linewidth=2)
plt.scatter(s3,n_max_mas(masa, s3), color= "red", linewidth=4)
plt.scatter(s6,n_max_mas(masa, s6), color= "red", linewidth=4)
plt.scatter(s9,n_max_mas(masa, s9), color= "red", linewidth=4, label="Geometric Model ")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig("f_nmax_elipsoid")
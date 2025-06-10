import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
import csv

archivos = ["4a7e.pdb", "5yce.pdb","1gwo.pdb","5u4y.pdb","3v03.pdb", "1gx4.pdb","1lfg.pdb","3qyt.pdb","1dgf.pdb","8cpv.pdb","5d4k.pdb", "HSA_autopsf.pdb"]

radio_gran = []
radio_mid  = []
radio_peq  = []
rog = []

for i in archivos:
    universe = mda.Universe(i)
    protein = universe.select_atoms("protein")
    masa = protein.total_mass()/1000 #kDa
    print("masa ", masa)
    rogg = protein.radius_of_gyration()/10
    rog.append(rogg)
    I_tensor = protein.moment_of_inertia()
    I_values = np.linalg.eigvals(I_tensor) /100 /1000 # Tres valores correspondientes a los ejes principales EN NANOMETROS^2 * kDa
    I_values.sort()  #ordena de menor a mayor
    radio1 = np.sqrt(I_values[0]*5/(2*masa) )           # el radio peuqeño, repetido b=c
    radio2 = np.sqrt((I_values[1]*5/masa)-radio1**2)    # el radio grande
    radio0 = np.sqrt((I_values[2]*5/masa)-radio1**2)    # el radio grande pero otra opción
    radios = []
    radios.append(radio1)
    radios.append(radio2)
    radios.append(radio0)
    radios.sort()
    #r_1.append(radios[1])  # esto es b, radio "grande" del elipsodie
    #r_2.append(radios[0])   # esto es a, radio "pequeño" del elipsodie
    #r_3.append(radios[2])
    #r_1.append(radio2)  # esto es b, radio "grande" del elipsodie
    #r_2.append(radio1)   # esto es a, radio "pequeño" del elipsodie
    #r_3.append(radio0)
    radio_gran.append(radios[1])
    radio_mid.append(radios[0])
    radio_peq.append(radios[2])

print("radiogrande  ", radio_gran)
print("radiomid  ", radio_mid)
print("radiopequeno ", radio_peq)

r1=7.5 #nm  radio de la nanoparticula de oro
def n_max_elip (a,b):
    alpha = math.atan((b)/(a+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    return nmax

insulin_nmax= n_max_elip(radio_peq[0], radio_gran[0])
myoglobin_nmax= n_max_elip(radio_peq[1], radio_gran[1])
horseradish_peroxidase_nmax= n_max_elip(radio_peq[2], radio_gran[2])
protein_A_nmax= n_max_elip(radio_peq[3], radio_gran[3])
BSA_nmax= n_max_elip(radio_peq[4], radio_gran[4])
galact_BSA_nmax = n_max_elip(radio_peq[5], radio_gran[5])
lactoferrin_nmax = n_max_elip(radio_peq[6], radio_gran[6])
transferrin_nmax = n_max_elip(radio_peq[7], radio_gran[7])
catalase_nmax = n_max_elip(radio_peq[8], radio_gran[8])
ferritin_nmax = n_max_elip(radio_peq[9], radio_gran[9])
polymeric_IgA_nmax = n_max_elip(radio_peq[10], radio_gran[10])
HSA_nmax = n_max_elip(radio_peq[11], radio_gran[11])





#datos del paper de 1987
insulin_weight = 5.88  # kDa
insulin_molecules = 200
myoglobin_weight = 17.00  # kDa
myoglobin_molecules = 87
horseradish_peroxidase_weight = 40.00  # kDa
horseradish_peroxidase_molecules = 61
protein_A_weight = 41.00  # kDa
protein_A_molecules = 60
BSA_weight = 68.00  # kDa
BSA_molecules = 39
galact_BSA_weight = 75.00  # kDa
galact_BSA_molecules = 15
lactoferrin_weight = 80.00  # kDa
lactoferrin_molecules = 35
transferrin_weight = 90.00
transferrin_weight = 69.99  # kDa
transferrin_molecules = 31
catalase_weight = 250.00  # kDa
catalase_molecules = 18
LDL_weight = 350.00  # kDa
LDL_molecules = 9
ferritin_weight = 450.00  # kDa
ferritin_molecules = 8
polymeric_IgA_weight = 500.00  # kDa
polymeric_IgA_molecules = 4
HSA_weight = 133.00  # kDa
HSA_molecules = 24

#funciones para calcular el n max
def radio(m):
    m1 = m*1000
    r2 = 0.066 * (m1**(1/3))
    return r2

def n_max(r2, r1):
    alpha = math.atan(r2/(r2+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    #print(nmax)
    return nmax

modelo =[]

radi = np.linspace(1,5.3)
radio_np = 7.5

#linea negra, modelo a partir de la masa
for i in radi:
    modelo.append(n_max(i, radio_np))
plt.figure(figsize=(8, 6))
plt.plot(radi, modelo, color="black", linewidth=1, label = "Original Model")



#pasar el radio a masa con la funcion radio
insulin_radius = radio(insulin_weight)
myoglobin_radius = radio(myoglobin_weight)
horseradish_peroxidase_radius = radio(horseradish_peroxidase_weight)
protein_A_radius = radio(protein_A_weight)
BSA_radius = radio(BSA_weight)
galact_BSA_radius = radio(galact_BSA_weight)
lactoferrin_radius = radio(lactoferrin_weight)
transferrin_radius = radio(transferrin_weight)
catalase_radius = radio(catalase_weight)
LDL_radius = radio(LDL_weight)
ferritin_radius = radio(ferritin_weight)
polymeric_IgA_radius = radio(polymeric_IgA_weight)
HSA_radius = radio(HSA_weight)

#prueba
insulin_radius = rog[0]
myoglobin_radius = rog[1]
horseradish_peroxidase_radius = rog[2]
protein_A_radius = rog[3]
BSA_radius = rog[4]
galact_BSA_radius = rog[5]
lactoferrin_radius = rog[6]
transferrin_radius = rog[7]
catalase_radius = rog[8]
#LDL_radius = rog[]
ferritin_radius = rog[9]
polymeric_IgA_radius = rog[10]
HSA_radius = rog[11]

#linea verde, modelo elipsoide
x_sim = np.array([
    insulin_radius, myoglobin_radius, horseradish_peroxidase_radius,protein_A_radius, BSA_radius, 
    galact_BSA_radius, lactoferrin_radius,transferrin_radius, catalase_radius, ferritin_radius, polymeric_IgA_radius, HSA_radius
]) 
y_sim = np.array([
    insulin_nmax, myoglobin_nmax, horseradish_peroxidase_nmax, protein_A_nmax, BSA_nmax, galact_BSA_nmax, lactoferrin_nmax, transferrin_nmax, catalase_nmax, ferritin_nmax, polymeric_IgA_nmax, HSA_nmax
])
coeficientes = np.polyfit(x_sim, y_sim, 2)  # Ajuste de curva
polinomio = np.poly1d(coeficientes)  # Crear la función del polinomio
curva_sim = polinomio(radi) 
plt.plot(radi, curva_sim, color="green", linewidth=1)

plt.scatter(insulin_radius, insulin_molecules, color='red')
plt.scatter(insulin_radius, insulin_nmax, color='green')
plt.annotate("Insulin", (insulin_radius, insulin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(myoglobin_radius, myoglobin_molecules, color='red')
plt.scatter(myoglobin_radius, myoglobin_nmax, color='green')
plt.annotate("Myoglobin", (myoglobin_radius, myoglobin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(horseradish_peroxidase_radius, horseradish_peroxidase_molecules, color='red')
plt.scatter(horseradish_peroxidase_radius, horseradish_peroxidase_nmax, color='green')
plt.annotate("Horseradish Peroxidase", (horseradish_peroxidase_radius, horseradish_peroxidase_molecules), fontsize=9, xytext=(0,5), textcoords='offset points')

plt.scatter(protein_A_radius, protein_A_molecules, color='red')
plt.scatter(protein_A_radius, protein_A_nmax, color='green')
plt.annotate("Protein A", (protein_A_radius, protein_A_molecules), fontsize=9, xytext=(-5,-8), textcoords='offset points')

plt.scatter(BSA_radius, BSA_molecules, color='red', label = "Experimental")
plt.scatter(BSA_radius, BSA_nmax, color='green', label = "Elipsoide")
plt.annotate("BSA", (BSA_radius, BSA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(galact_BSA_radius, galact_BSA_molecules, color='red')
plt.scatter(galact_BSA_radius, galact_BSA_nmax, color='green')
plt.annotate("Galact. BSA", (galact_BSA_radius, galact_BSA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(lactoferrin_radius, lactoferrin_molecules, color='red')
plt.scatter(lactoferrin_radius, lactoferrin_nmax, color='green')
plt.annotate("Lactoferrin", (lactoferrin_radius, lactoferrin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(transferrin_radius, transferrin_molecules, color='red')
plt.scatter(transferrin_radius, transferrin_nmax, color='green')
plt.annotate("Transferrin", (transferrin_radius, transferrin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(catalase_radius, catalase_molecules, color='red')
plt.scatter(catalase_radius, catalase_nmax, color='green')
plt.annotate("Catalase", (catalase_radius, catalase_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(LDL_radius, LDL_molecules, color='red')
plt.annotate("LDL", (LDL_radius, LDL_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(ferritin_radius, ferritin_molecules, color='red')
plt.scatter(ferritin_radius, ferritin_nmax, color='green')
plt.annotate("Ferritin", (ferritin_radius, ferritin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(polymeric_IgA_radius, polymeric_IgA_molecules, color='red')
plt.scatter(polymeric_IgA_radius, polymeric_IgA_nmax, color='green')
plt.annotate("Polymeric IgA", (polymeric_IgA_radius, polymeric_IgA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(HSA_radius, HSA_molecules, color='red')
plt.scatter(HSA_radius, HSA_nmax, color='green')
plt.annotate("HSA", (HSA_radius, HSA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')


plt.xlabel("Protein Radius (nm)")
plt.ylabel("Molecules per nanoparticle")
plt.title("Protein Binding ")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.legend()
plt.savefig("g_grafica_comparacion")

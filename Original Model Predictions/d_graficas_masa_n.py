import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
import csv


# m es la masa de la proteina en kDat
# r1 es el radio de la nanoparticula
def n_max (m, r1):
    #print(m)
    #m1 = (m / (1.66054 * 10**(-21)) ) #* (6.02214076*10**23)
    #print(m1)
    m1 = m*1000
    r2 = 0.066 * (m1**(1/3))
    #print(r2)
    alpha = math.asin(r2/(r2+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    #print(nmax)
    return nmax


modelo =[]

masa = np.linspace(5,500)
radio_np = 7.5

for i in masa:
    modelo.append(n_max(i, radio_np))



plt.figure(figsize=(8, 6))
#plt.xlim(230, 530)
#plt.ylim(0, 35) 
plt.plot(masa, modelo, color="black", linewidth=1, label="Geometrical Model")

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
transferrin_weight = 90.00  # kDa
transferrin_molecules = 31
catalase_weight = 250.00  # kDa
catalase_molecules = 18
LDL_weight = 350.00  # kDa
LDL_molecules = 9
ferritin_weight = 450.00  # kDa
ferritin_molecules = 8
polymeric_IgA_weight = 500.00  # kDa
polymeric_IgA_molecules = 4




# Graficar cada punto con su etiqueta al lado
plt.scatter(insulin_weight, insulin_molecules, color='red')
plt.annotate("Insulin", (insulin_weight, insulin_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(myoglobin_weight, myoglobin_molecules, color='red')
plt.annotate("Myoglobin", (myoglobin_weight, myoglobin_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(horseradish_peroxidase_weight, horseradish_peroxidase_molecules, color='red')
plt.annotate("Horseradish Peroxidase", (horseradish_peroxidase_weight, horseradish_peroxidase_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(protein_A_weight, protein_A_molecules, color='red')
plt.annotate("Protein A", (protein_A_weight, protein_A_molecules), fontsize=13, xytext=(-40,-10), textcoords='offset points')

plt.scatter(BSA_weight, BSA_molecules, color='red')
plt.annotate("BSA", (BSA_weight, BSA_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(galact_BSA_weight, galact_BSA_molecules, color='red')
plt.annotate("Galact. BSA", (galact_BSA_weight, galact_BSA_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(lactoferrin_weight, lactoferrin_molecules, color='red')
plt.annotate("Lactoferrin", (lactoferrin_weight, lactoferrin_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(transferrin_weight, transferrin_molecules, color='red')
plt.annotate("Transferrin", (transferrin_weight, transferrin_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(catalase_weight, catalase_molecules, color='red')
plt.annotate("Catalase", (catalase_weight, catalase_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(LDL_weight, LDL_molecules, color='red')
plt.annotate("LDL", (LDL_weight, LDL_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(ferritin_weight, ferritin_molecules, color='red')
plt.annotate("Ferritin", (ferritin_weight, ferritin_molecules), fontsize=13, xytext=(5,5), textcoords='offset points')

plt.scatter(polymeric_IgA_weight, polymeric_IgA_molecules, color='red', label= "Experimental measurements")
plt.annotate("Polymeric IgA", (polymeric_IgA_weight, polymeric_IgA_molecules), fontsize=13, xytext=(-70,-10), textcoords='offset points')


plt.xlabel(r"M$_r$ (kDa)", fontsize=14)
plt.ylabel(r"N$_{max}$ of Proteins in Gold Nanoparticles of r= 7.5 nm", fontsize=14)
plt.ylabel(r"N$_{max}$", fontsize=14)
#plt.title("Proof of M. Soloviev et al (2022) model")

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.legend(fontsize=14)
plt.savefig("g_graficas_masa")


import MDAnalysis as mda
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
import csv

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


def radio(m):
    m1 = m*1000
    r2 = 0.066 * (m1**(1/3))
    return r2

def n_max(r2, r1):
    alpha = math.asin(r2/(r2+r1))
    nmax = math.floor((2 * math.sqrt(3)* math.pi) / (3 * (alpha ** 2)))
    #print(nmax)
    return nmax

modelo =[]

radi = np.linspace(1,5.5)
radio_np = 7.5

for i in radi:
    modelo.append(n_max(i, radio_np))


plt.figure(figsize=(8, 6))
plt.plot(radi, modelo, color="black", linewidth=1)


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

# Graficar cada punto con su etiqueta al lado usando radius en lugar de weight
plt.scatter(insulin_radius, insulin_molecules, color='red')
plt.annotate("Insulin", (insulin_radius, insulin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(myoglobin_radius, myoglobin_molecules, color='red')
plt.annotate("Myoglobin", (myoglobin_radius, myoglobin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(horseradish_peroxidase_radius, horseradish_peroxidase_molecules, color='red')
plt.annotate("Horseradish Peroxidase", (horseradish_peroxidase_radius, horseradish_peroxidase_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(protein_A_radius, protein_A_molecules, color='red')
plt.annotate("Protein A", (protein_A_radius, protein_A_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(BSA_radius, BSA_molecules, color='red')
plt.annotate("BSA", (BSA_radius, BSA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(galact_BSA_radius, galact_BSA_molecules, color='red')
plt.annotate("Galact. BSA", (galact_BSA_radius, galact_BSA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(lactoferrin_radius, lactoferrin_molecules, color='red')
plt.annotate("Lactoferrin", (lactoferrin_radius, lactoferrin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(transferrin_radius, transferrin_molecules, color='red')
plt.annotate("Transferrin", (transferrin_radius, transferrin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(catalase_radius, catalase_molecules, color='red')
plt.annotate("Catalase", (catalase_radius, catalase_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(LDL_radius, LDL_molecules, color='red')
plt.annotate("LDL", (LDL_radius, LDL_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(ferritin_radius, ferritin_molecules, color='red')
plt.annotate("Ferritin", (ferritin_radius, ferritin_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.scatter(polymeric_IgA_radius, polymeric_IgA_molecules, color='red')
plt.annotate("Polymeric IgA", (polymeric_IgA_radius, polymeric_IgA_molecules), fontsize=9, xytext=(5,5), textcoords='offset points')

plt.xlabel("Protein Radius (nm)")
plt.ylabel("Molecules per Particle at Saturation")
plt.title("Protein Binding ")

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("g_graficas_radio")

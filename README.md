# Repository TFG -  Mario Martínez Angulo

This repository contains the Python codes used and the molecular dynamics simulations created during my undergraduate thesis, as well as the results, graphs, and images I obtained. The codes are separated into different files according to the purpose and section of each one.

- Radius of gyration. Here you can find the codes for calculating the radius of gyration of HSA and the calculation of the maximum number of proteins on the surface through the method described by [M. Soloviev et al., 2021](https://doi.org/10.1016/j.jcis.2021.07.072). You can also find images made by VMD of the HSA protein and the graphs of the code.

- Maximum Distance: Codes that perform the calculation for the maximum distance between two atoms in each frame of the simulation. Using the same geometric method as before, we calculate the maximum number of proteins on the nanoparticle (NP), but now considering this distance. Still considering a espheric protein with this maximun distances as the diameter.

- Conformational Changes: Another approach involves examining how the protein undergoes conformational changes when it adsorbs onto the surface of a nanoparticle, and how this movement affects the maximum number of proteins.

- Inertia Radius: Calculations of the three principal moments of inertia and their corresponding radius, approximating the protein as a sphere.

- Ellipsoid. The geometric model and the calculation of the moments of inertia are adjusted to an ellipsoidal shape for the protein.

- Original Model Predictions: Here you'll find graphical recreations of the original results from [M. Soloviev et al., 2021](https://doi.org/10.1016/j.jcis.2021.07.072) 

- Comparation: comparative analysis between their model and the ellipsoidal approximation developed in this work.

- BUBBLES: Configuration scripts used for setting up BUBBLES simulations, including all relevant parameter files. Note that the core simulation code remains unchanged from the official BUBBLES release; only input configurations were adapted for this study.

- Protein Structures and HSA: This folder contains the .pdb files of the proteins used, including the HSA structures derived from molecular dynamics simulations.

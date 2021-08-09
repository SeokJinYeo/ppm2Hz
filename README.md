# ppm2Hz : Fast and accurate Matlab based 'susceptibility-induced B0 inhomogeneity' calculation programs
ppm2Hzs are based on generalized susceptibility voxel convolution(gSVC) method which is rapid and artifact-free[1].
For application in various cases, the method was extended to arbitrary orientation and spatially varying applied field case[2].
Additionally another static magnetic field perturbation calculation method, k-space-discretized(KD) is coded(ppm2HzKD). If you want more details of method, see the references below.

# Usage examples
- example of ppm2Hz applied at ellipsoid model with arbitrary direction 
- example of spatially varying applied field

# Versions
- ppm2Hz : Conventional cases when direction of object's magnetization(n) and the direction in which its magnetic field(l) are coincide to z-direction(n=l=z).
- nppm2Hz : When two directions are coincede to arbitrary direction(n=l).
- ppm2Hznl : When two directions are set arbitrary.
- ppm2Hznl2 : When applied field is spatially varying and direction of perturbation is arbitrary.

# Common input parameters
- chi : source susceptibility model 
- dr : voxel size [m]
- B0 : strength of applied magnetic field [T]
- n : direction of applied magnetic field(magnetization)
- l : direction of perturbation field (n and l are unit vectors)

# gSVC input parameters
- r : distance between source and target grid [m]
- t : target grid size

# KD input parameters
- ff : fine-grain factor
- bf : zero-filling factor


# References
This work has been published in MRM.
[1] Seung‐Kyun Lee, Seon‐Ha Hwang, Ji‐Seong Barg, Seok‐Jin Yeo. (2018). “Rapid, theoretically artifact‐free calculation of static magnetic field induced by voxelated susceptibility distribution in an arbitrary volume of interest”, Magnetic resonance in medicine.
[2] Seok-Jin Yeo, So-Hee Lee, Jang-Yeon Park, Seung-Kyun Lee. (2021). “Susceptibility voxel convolution for calculation of MRI static magnetic field perturbation from uniformly magnetized objects in arbitrary orientations.”, Magnetic resonance in medicine(now revisioning).




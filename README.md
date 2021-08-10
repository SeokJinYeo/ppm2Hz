# ppm2Hz : Fast and accurate Matlab based 'susceptibility-induced B0 inhomogeneity' calculation programs
ppm2Hzs are based on generalized susceptibility voxel convolution(gSVC) method which is rapid and artifact-free[1].
For application in various cases, the method was extended to arbitrary orientation and spatially varying applied field case[2].
Additionally another static magnetic field perturbation calculation method, k-space-discretized(KD) is coded(ppm2Hz_KD). If you want more details(theory, applications) of method, see the references below.

# Usage examples
- [examples][elliplink] of ppm2Hz applied at ellipsoid model

[elliplink]: https://github.com/SeokJinYeo/ppm2Hz/tree/main/Example "example"

# Versions
Requiring input parameters are varying in each version, users need to check them at each code.
- ppm2Hz : Conventional cases when direction of object's magnetization(n) and the direction in which its magnetic field(l) are coincide to z-direction(n=l=z).
- nppm2Hz : When two directions are coincede to arbitrary direction(n=l).
- ppm2Hznl : When two directions are set arbitrary.
- ppm2Hznl2 : When the applied field is spatially varying and direction of perturbation is arbitrary. We need scalar field of each component of applied vector magnetic field(Bx,By,Bz) .Those scalar fields need to have the same grid and voxel size with input susceptibility model.

# Common input parameters
- chi : source susceptibility model 
- dr : voxel size [m]
- B0 : strength of applied magnetic field [T]
- n : direction of applied magnetic field(magnetization)
- l : direction of perturbation field (n and l are unit vectors)

# gSVC input parameters
- r : distance between source and target grid [m]
- t : target grid size
- **ex) b0 = ppm2Hznl(chi,dr,r,t,n,l,B0);**

# KD input parameters
- ff : fine-grain factor
- bf : zero-filling factor
- **ex) b0 = ppm2Hznl_KD(chi,dr,ff,bf,n,l,B0);**

# References
This work has been published in MRM.
[1] Seung‐Kyun Lee, Seon‐Ha Hwang, Ji‐Seong Barg, Seok‐Jin Yeo. (2018). “Rapid, theoretically artifact‐free calculation of static magnetic field induced by voxelated susceptibility distribution in an arbitrary volume of interest”, Magnetic resonance in medicine.
[1] Seok-Jin Yeo, So-Hee Lee, Seung-Kyun Lee. (2021). “Rapid calculation of static magnetic field perturbation generated by magnetized objects in arbitrary direction”. Magnetic resonance in medicine(예정). Page.




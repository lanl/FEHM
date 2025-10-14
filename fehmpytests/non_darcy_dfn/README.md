---
title : fehmpytests non_darcy_dfn
---

# non_darcy_dfn

**Test Darcy and Non-Darcy Water into a Single dfn fracture**


These are a set of tests comparing Darcy and non-Darcy flow simulations using FEHM on a single fracture, modeled as a 2D rectangle with a thickness defined as the aperture. The input files for these simulations were created using dfnWorks, although dfnWorks is not required for running these FEHM tests. 

This test consists of 9 representative cases selected from the complete set of consisting of 441 simulation runs, which include:
    7 different beta values for the ndar macro (ranging from 1e-10 to 1e-6)
    3 different apertures (1e-3, 5e-3, and 1e-4)
    21 different pressure differential (dp) values, where dp1e01 indicates a pressure change of 10 MPa.

Python scripts were used to create run directories and files for FEHM runs. The FEHM input decks for each run were modified to control the darcy or non-darcy run. The ndar macro defines the beta value used, changes in pressure differential (dp) are made by modifying the pressures in the flow macro zones. The aperture is adjusted bynaming one of the three different area/volume stor file in the fehmn.files list.

The graphs illustrate water discharge (Q) plotted against the change in pressure (dp) as a function of aperture. It is essential to note that Q is influenced by factors such as aperture size, permeability, and cross-sectional area, rather than just velocity/Re. The Q values represented in the plots are derived from the last time-step recorded in the .out file and are converted to cubic meters per second (m³/s) by dividing by the water density of 998 kg/m³. This density, denoted as rho, is defined in the FEHM input deck, assuming a constant temperature of 20°C. 


<p> 
 <a href="_information/cplot_aperture_1.0mm.png"> <img width="500" src="_information/cplot_aperture_1.0mm.png"> </a> 
</p>

Run 9 tests with the following:

python fehmpytests.py ../src/xfehm non_darcy_dfn 

A1B2_dp1e00 A1B5_dp1e00 A1Darcy_dp1e00
A2B2_dp1e00 A2B5_dp1e00 A2Darcy_dp1e00
A3B2_dp1e00 A3B5_dp1e00 A3Darcy_dp1e00

A1 = aperture e1e-3 = 0.001
A2 = aperture e5e-4 = 0.0005
A3 = aperture e1e-4 = 0.0001
B2 = ndar beta = 1.0e+02
B5 = ndar beta = 1.0e+05

LANL server locations:
nondarcy runs: /scratch/er/mdalal/nondarcy/batch_runs
darcy runs: /scratch/er/mdalal/nondarcy/darcy_batch_runs


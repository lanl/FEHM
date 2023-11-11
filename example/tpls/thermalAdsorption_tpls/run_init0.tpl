ptf %
title: REGOLITH ADSORPTION. (Establish air-static and geothermal gradient).
sol
1 -1
#----------NGAS/COND-----------
# Need ngas/cond to assign custom EOS
#1 0 0   1  0.5 0.5 0 
#1 0 0   1  %thermal_cond% %thermal_cond% 0 
cond
1 0 0    0.5  0.5  0.5

#----------EOS-----------------
eos
file
../macros/eos_mars_co2.dat
#----------INIT----------------
# For fully saturated water flow problems (otherwise use pres)
# Surface temp (upshifted for FEHM), with linear geothermal gradient 
# for depths between 0 and 500 m
#--- Also requires: cond
# PEIN    TIN  TIN1  GRADS1  DEPTH  TIN2  GRAD2  QUAD
init
%P_initial%  0.0  %T_0%  0.  1.e4   %T_0%  0.   0. 
# 0.70618  0.0  %T_0%  0.  1.e4   %T_0%  0.   0. 
# 0.070618  0.0  %T_0%  0.  1.e4   %T_0%  0.   0. 
# 0.070618  0.0  %T_0%  %geothermal_gradient%  1.e4   %T_0%  %geothermal_gradient%   0. 
# 0.06  0.0  0.05  0. 0.  0.05  0. 0. 
#----------ZONES---------------
# 100 - top
#   9 - fracture
#  12 - obs nodes at top and bottom
# zone
# file
# ../macros/top.zone
zone
file
../macros/obs_topbot.zone
#----------NODE----------------
# node
# 2  
# -88 -89
# 0. 0. 0. 
# 0. -50. 0.
# node
# block
# -9 0 0 
node
block
-12 0 0 
 
#---------------OUTPUT----------------------
hist
days 1.e20 10000 
pressure
density
temperature
end
# cont
# tec 1e20  1.e8 days
# xyz
# material
# liquid
# pressure
# endcont
# # geo
#---------------Time/Ctrl-------------------
# NOTE: convergence/out-of-bounds issues when G3 and TMCH tol
time
 1.000 365.d8  10000000  50  2007 1

ctrl
30 1.e-05 30 100 gmres
1 0 0 2

1.00000 2 1.0
10 1.6 1e-9  365.d6
1 -1
iter
1e-5  1e-5 1e-2  -1e-3  1.2
0       0       0       5       14000.
#----------------PERM/ROCK/RLP---------------
zone
file
%matzone_file%
perm
file
../macros/perm.dat
rock
file
../macros/rock.dat
#----------------BOUN------------------------
zone
file
../macros/top.zone
# zonn
# file
# ../macros/bottom.zonn
boun
file 
../macros/isothermal.boun
# ../macros/geo-grad_air-static.boun
stop


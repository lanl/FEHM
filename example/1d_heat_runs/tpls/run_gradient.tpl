ptf %
title: Mars heatflow init0 run. 
sol
1 -1
#----------NGAS/COND-----------
# Need ngas/cond to assign custom EOS
cond
1 0 0   1  2.5 2.5 0 

#----------EOS-----------------
eos
file
../macros/eos_mars_co2_test.dat
#----------INIT----------------
# For fully saturated water flow problems (otherwise use pres)
# Surface temp = 0.0°C (FEHM), with linear geothermal gradient = 0.0045 °C/m
# for depths between 0 and 500 m
#--- Also requires: cond
# PEIN    TIN  TIN1  GRADS1  DEPTH  TIN2  GRAD2  QUAD
init
0.070618  0.0  %T_0%  %geothermal_gradient%  1.e4   %T_0%  %geothermal_gradient%   0. 
# 0.06  0.0  0.05  0. 0.  0.05  0. 0. 
#----------ZONES---------------
# 100 - top
#   9 - fracture
# zone
# file
# ../macros/top.zone
#----------NODE----------------
# node
# 2  
# -88 -89
# 0. 0. 0. 
# 0. -50. 0.
# node
# block
# -9 0 0 
# 
node
file
../macros/nodeline.dat
#---------------OUTPUT----------------------
hist
days 1.e20 10000 
pressure
temperature
end
# density
# cont
# tec 1e20  1.e8 days
# xyz
# material
# liquid
# pressure
# endcont
# # geo
#---------------Time/Ctrl-------------------
time
 1.000 365.d8  10000000  50  2007 1

ctrl
30 1.e-05 30 100 gmres
1 0 0 2

1.00000 2 1.0
10 1.6 1e-9  365.d7
1 0 
iter
1e-5  1e-5 1e-4  -1e-6  1.2
0       0       0       5       14000.
#----------------PERM/ROCK/RLP---------------
zone
file
../macros/fracture.zone
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
zonn
file
../macros/bottom.zonn
boun
file 
../macros/air_static.boun
# boun
# model 1
# ti
# 2
# 0
# 1e10
# pw
# 0.070618
# 0.070618
# t
# 0.00
# 0.00
# model 2
# ti
# 2
# 0.0  1e20
# t
# %T_bot%  %T_bot%
# 
# -100 0 0  1
# -200 0 0  2
# 
# #----------------FLOW------------------------
# flow
# -100 0 0   1e-8 -0.05 1e3
# 1 1 1     -1e-8 -0.05 1e3   
# 
stop

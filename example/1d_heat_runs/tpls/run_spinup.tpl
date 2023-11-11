ptf %
title: Mars periodic heat flow spinup run. 
# sol: NTT INTG
# If NTT is 1,  it's a coupled solution
# If NTT is -1, it's a heat transfer only solution
sol
-1 -1
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
0.070618  0.0  %T_0%  0.0  1.e4   %T_0%  0.0   0. 
# 0.070618  0.0  %T_0%  %geothermal_gradient%  1.e4   %T_0%  %geothermal_gradient%   0. 
# 0.06  0.0  0.05  0. 0.  0.05  0. 0. 
#----------ZONES---------------
# 100 - top
#   9 - fracture
zone
file
../macros/top.zone
#----------NODE----------------
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
#---------------Time/Ctrl-------------------
time
 1e-9 %sim_time%  10000000   1 1 1 0. 

ctrl
30 1.e-07 30 100 gmres
1 0 0 2

1.00000 2 1.0
10 1.6 1e-9  1e-1
1 0 
iter
1e-5  1e-5 1e-4  -1e-6  1.2
0       0       0       5       14000.
#----------------PERM/ROCK/RLP---------------
# zone
# file
# ../macros/fracture.zone
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
boun
file 
%boun_file_no_geo%
stop


ptf %
title: REGOLITH ADSORPTION. Mars baro pumping w/ adsorption. 
sol
1 -1
#----------NGAS/COND-----------
# Need ngas/cond to assign custom EOS
# 1 0 0   1  0.5 0.5 0 
# 1 0 0   1  thermal_cond thermal_cond 0 
cond
1 0 0    0.5 0.5 0.5

#----------EOS-----------------
eos
file
../macros/eos_mars_co2.dat
#----------ZONES---------------
# 100 - top
#   9 - fractures
zone
file
../macros/top.zone
zone
file
%matzone_file%
#----------NODE----------------
#will need to fix this
node
file
nodeline.dat
#---------------OUTPUT----------------------
hist
days 1.e20 10000 
pressure water
density water
viscosity water
temperature
concentration
velocity
zflux
global mass
end
#----CONT (for occasional run.fin dumps)----
# Dump a run.fin file every 1000 days
cont
tec 1e30 1000. 
material
concentration
endcont
#---------------cflx------------------------
# output solute fluxes/sources/sinks for top of domain (only does aqueous methane)
zone
file
../macros/top.zone
cflx
1
100
#---------------Time/Ctrl-------------------
# 2712.41 sols = 2786.96 days
#0.0000001 %sim_time%  100000000  1 1 1 0. 
time
0.0000001 %sim_time%  100000000  1 1 1 0. 

ctrl
30 1.e-05 30 100 gmres
1 0 0 2

1.00000 2 1.0
10 1.6 1e-9  1e-1
1  0 
iter
1e-5  1e-5 1e-2  -1e-3  1.2
0       0       0       5       1.e11 
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
# 100 - top nodes 
# 200 - bottom nodes
# 8 - contaminated zone
zone
file
../macros/top.zone
zonn
file
../macros/contaminated_layer.zonn
# zonn
# file
# ../macros/bottom.zonn
boun
file
../macros/%boun_file%
#------------TRAC----------------------------
# 8 - contaminated zone
zone
file
../macros/contaminated_layer.zone
trac
file
../macros/trac.dat
#-------------RXN-------------------------
# 100 - top zone
# 11  - constant regolith (no temp change in sorption)
zone
file
../macros/top.zone
zonn
file
../macros/constant_regolith.zonn
# Temperature-dependent Langmuir adsorption
rxn
file
../macros/rxn.dat
# ../macros/rxn_file.dat
stop

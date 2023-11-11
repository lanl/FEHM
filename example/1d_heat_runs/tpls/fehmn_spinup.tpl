ptf %
root: run
input: run.dat
grida: %mesh_file% 
rsto: run.fin
error: run.err
check: run.chk
outp: run.out

none
0

#Don't restart (don't want geotherm gradient for analytical soln)
rsti: %init0_dir%/run.fin
# grida: /project/gas_seepage/jportiz/mars/mesh/2d_planar/fracture_5x50m/grid.inp
# storo: stor_file


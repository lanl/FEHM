---
title : heat2d_quad
layout : page_getting-started
hero_height: is-hidden
---

# heat2d_quad

This is heat2d_quad in the VV Test Suite, output files are written to heat_flux directory. 
Files have been modified to read and write files under the heat2d_quad directory and there is no heat_flux.

From readme file in VV test directory:
<pre>
Checked with the results in heat2d_quad.00003_hf_node.dat (dated Sept 23, 2014),
it agrees with analytical solution as calculated below:

Heat copnduction:
q = k*delta-T/delta-L

K=2.7 W/m/degK
deltaT = 100 deg C
delta-L = 0.5m
q = 2.7*100/.5 =5.4e+2 (W/m/degK)*(degC)/m = 5.4e+2 W/m^2
= 5.4e-4 MW/m^2
= 0.54e-3 MW/m^2
</pre>


Test Directory: [FEHM/fehmpytests/heat2d_quad](https://github.com/lanl/FEHM/tree/master/fehmpytests/heat2d_quad)


### Example File heat2d_quad.in 
<pre>
***** 2-D Heat Conduction Model *****
# Added following to be consistent with original test
nfinv
fhot
node
  11
  1 12 56 100 111
  61
  11 22 66  110 121
cont
tecplot 1000 10.
geom
temperature
heatflux
end
sol
  -1    -1
init
  10.    0.    200.   0.   0.   200.   0.   0.
zone
  1
  0.       0.500000       0.500000       0.
  0.       0.             0.500000       0.500000

rock
   -1    0    0    2700.     1000.     0.      1.0

cond
   -1    0    0    2.7e-00   2.7e-00   2.7e-00

perm
   -1    0    0    1.e-30    1.e-30    1.e-30  0.  0.  0.

zone
  2
  -0.001 +0.001 +0.001 -0.001
  -0.001 -0.001 5.001 5.001
  3
  0.499 0.501 0.501 0.499
  -0.001 -0.001 5.001 5.001

flow
  -2 0 0     10.00   -100.00  1.e03
  -3 0 0     10.00   -200.00  1.e03

time
  0.005   4.00   1000   10   1989   04
  .25 .005 1. 10

ctrl
  40   1.e-04   08
  -1    0    0    1

  1.0   0.0   1.0
  10   1.0   0.00005 0.005
  1   0
stop
</pre>

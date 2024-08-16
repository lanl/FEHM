---
title : wvtest
layout : page_getting-started
hero_height: is-hidden
---

# wvtest

**Evaporation Heatpipe Dry-out Test**

This is a test of the ability of the code to do dry-out by evaporation.

The problem is 1 m long with 5 elements and 12 nodes,
the left side is fixed at 35 C and the right side is fixed at 10 C.

Initial saturation is 0.25  (left) to 0.15 (right) in steps of 0.02

porosity = 0.8   tortuosity(adif) = 0.66

Capillary forces are turned off so that as the water is moved to
the cold end there is no way for it to flow back to the hot end.
Basically this is a one-way heat pipe.


Test Directory: [FEHM/fehmpytests/wvtest](https://github.com/lanl/FEHM/tree/master/fehmpytests/wvtest)


### Example File box.dat 
<pre>
#   Test of Water Vapor Diffusion Water Vap Tracer
#
#  1 Vapor Tracer   molecular diffusion  set to zero!
#
#   driven  by delta T= 25C
#   No capillary forcing
#
cont
tecplot 30000  1000.
geom
material
liquid
conc
veloc
vapor
temp
pres
sat
endavs
hist tecplot
mpa
saturation
end
sol
  1    -1
rock
 1   0 0   1000.  1000.  0.8
0
pres
  1   2      1      0.10     0.25   2
  3   4      1      0.10     0.23   2
  5   6      1      0.10     0.21   2
  7   8      1      0.10     0.19   2
  9  10      1      0.10     0.17   2
 11  12      1      0.10     0.15   2
  0
adif
0.66
ngas
3
 1  2   1       -35.0
 3  4   1       -30.0
 5  6   1       -25.0
 7  8   1       -20.0
 9 10   1       -15.0
11 12   1       -10.0

1       0       0       0.     0.

1 0 0   0.00  0.

rlp
 2  0. 0. 0. 0.

 1 0 0 1

node
6
 1 3 5 7 9 11
perm
   1    0 0    1.e-12   1.e-12  1.e-12
   0
hflx
 1   2   1    35.  1e6
11  12   1    10.  1e6

flxo
5
1 3
3 5
5 7
7 9
9 11
finv
cond
 1   0 0  10.0  10.0    10.0
0
time
  0.001  5000.  55000    001  0000   01
  0.0    0.0    0.0    0.0
ctrl
  -8   1.e-07  8
  1  0 0     4
    0    0    0    0
  1.0   0   1.0
  07   1.4   1.e-7  5.
  1  0
trac
0.0     1.0      1.e-7   1.0
0       1.0e6    1.0e6   0
50      1.4      1.e-3  5.0
1
ldsp
-1
0  0 0 1  2 2.0e-5 0. 0.

1 0 0 1

1 0 0       1.


stop

</pre>

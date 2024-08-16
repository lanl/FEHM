---
title : heat3d
layout : page_getting-started
hero_height: is-hidden
---

# heat3d

**test heat3d**

3-D Heat Conduction Problem with
Comparison of Model and Analytical Solution for Temperature vs Time

5 tests are run, see heat3d_mix.files  heat3d_quad.files  heat3d_ref.files  heat3d_tets.files  heat3d_tri.files

Test Directory: [FEHM/fehmpytests/heat3d](https://github.com/lanl/FEHM/tree/master/fehmpytests/heat3d)


### Example File heat3d.in 
<pre>
***** 3-D Heat Conduction Model *****
# To be consistent with earlier versions
nfinv
node
  1
  1321
cont
avs 1000 10.
geom
temperature
end
sol
  -1    -1
init
  10.    0.    200.   0.   0.   200.   0.   0.
zone
  1
          0.       0.500000       0.500000             0.
          0.       0.500000       0.500000             0.
          0.             0.       0.500000       0.500000
          0.             0.       0.500000       0.500000
          0.             0.             0.             0.
    0.500000       0.500000       0.500000       0.500000
  2
          0.       0.4800000       0.4800000             0.
          0.       0.4800000       0.4800000             0.
          0.             0.       0.4800000       0.4800000
          0.             0.       0.4800000       0.4800000
          0.             0.             0.             0.
    0.4800000       0.4800000       0.4800000       0.4800000

rock
   -1    0    0    2700.     1000.     0.      1.0
   -2    0    0    2700.     1000.     0.      1.0

cond
   -1    0    0    2.7e-00   2.7e-00   2.7e-00
   -2    0    0    2.7e-00   2.7e-00   2.7e-00

perm
   -1    0    0    1.e-30    1.e-30    1.e-30  0.  0.  0.
   -2    0    0    1.e-30    1.e-30    1.e-30  0.  0.  0.

flow
   -1    0    0     10.00   -100.00  1.e03

time
  0.005   3.00   1000   10   1989   04
  .25 .005 1. 10

ctrl
  40   1.e-04   08
   -1    0    0    1
   -2    0    0    1

  1.0   0.0   1.0
  10   1.0   0.00005 0.005
  0   0
stop
</pre>

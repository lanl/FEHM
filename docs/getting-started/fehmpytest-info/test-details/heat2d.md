---
title : heat2d
layout : page_getting-started
hero_height: is-hidden
---

# heat2d

**Test heat2d**

2-D Heat Conduction Problem with 3-node Triangles. Comparison of Model and Analytical Solution for Temperature vs Position and Comparison of Model and Analytical Solution for Temperature vs Time. 



Test Directory: [FEHM/fehmpytests/heat2d](https://github.com/lanl/FEHM/tree/master/fehmpytests/heat2d)


### Example File heat2d.in 
<pre>
***** 2-D Heat Conduction Model *****
# To be consistent with earlier versions
nfinv
node
  1
  111
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
  0.       0.500000       0.500000       0.
  0.       0.             0.500000       0.500000

rock
   -1    0    0    2700.     1000.     0.      1.0

cond
   -1    0    0    2.7e-00   2.7e-00   2.7e-00

perm
   -1    0    0    1.e-30    1.e-30    1.e-30  0.  0.  0.

flow
    1    11    1     10.00   -100.00  1.e03
   11   121   11     10.00   -100.00  1.e03

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

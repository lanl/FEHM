---
title : head
layout : page_getting-started
hero_height: is-hidden
---

# head

**Test Head Pressure Problem**

Comparison of Head and Pressure Formulation for Pressure vs Position at Time  365.000.
Compares the generated contour files to old contour files that are known to be correct. 

Test Directory: [FEHM/fehmpytests/head](https://github.com/lanl/FEHM/tree/master/fehmpytests/head)


### Example File head.dat
<pre>

Head versus Pressure Test - Head formulation
#***************************************************
cont 
avs     5000    1.00000e+19
geom
liquid 
head
pressure
formatted 
endavs
#***************************************************
head
airw
3
30.0 0.1
#bous
#1
#***************************************************
zone
00001 
0     0  100  100   0   0 100 100
100   0    0  100 100   0   0 100
100 100  100  100   0   0   0   0

#***************************************************
perm
   1 0 0 1e-13 1e-13 1e-14

pres
1 0 0 150. 1. 1

rock                              
  1    0     0     2.50E+03     1.00E+03    .1     

#******************************************************
zone
00060
 -1  -1   0   0  -1  -1   0   0
101  -1  -1 101 101  -1  -1 101
101 101 101 101  -1  -1  -1  -1
00061
100 100 101 101 100 100 101 101
101  -1  -1 101 101  -1  -1 101
101 101 101 101  -1  -1  -1  -1 
00072
 -1  -1 101 101   -1  -1 101 101 
101  -1  -1 101  101  -1  -1 101
101 101 101 101  100 100 100 100 

#******************************************************
flow
 -72 0 0 140. 1 1e6

# -60 0 0 140. 1 1e6
#***************************************************
time
1e-3 365 500000 3 1999 5 0.00000
 
#***************************************************
ctrl
 15 1.00000e-05 20
 1 0 0 2

 1.00000 3 1.00000
 10 2.00000 1.000e-08 1e20
 0 0
sol
 1 -1
iter
 1.00000e-05 1.00000e-05 0.00100000 -1.00000e-6 1.1
13       0       0       5       3600.
stop

</pre>

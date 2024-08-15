---
title : heat_pipe
layout : page_getting-started
hero_height: is-hidden
---

# heat_pipe

**Test the Heat Pipe Problem**

Compares the generated output files with the old files known to be correct. All values are tested.

RLPM subcases were excluded because they caused FEHM to produce too many negative volumes.

Test Directory: [FEHM/fehmpytests/heat_pipe](https://github.com/lanl/FEHM/tree/master/fehmpytests/heat_pipe)


### Example File heat_pipe_fdm.dat 

<pre>

title: 1-d heat pipe calculation                5-04-2001
#************************************************************************75
zone 
1 
0. 1. 1. 0.
0. 1. 1. 0.
0. 0. 1. 1.
0. 0. 1. 1.
0. 0. 0. 0.
-0.1 -0.1 -0.1 -0.1

node
10
1 2 3 10 20 30 40 48 49 50 
#************************************************************************75
grad
1
1 0. 3 2 10.  -150.
#************************************************************************75
perm  in m**2
1    0   0   3.9e-14 3.9e-14 3.9e-14      

#************************************************************************75
rlp
        3  0. 1.0  3.6    1.56   2.  0.0100 

        1      0       0       1

#************************************************************************75
rock  density,heatcapacity,porosity
        1      0       0       1500.   1100.   0.499

#************************************************************************75
cond Thermal conductivity
1      0       0     6.e-01   6.e-01   6.e-01

#************************************************************************75
pres       initial pres sat
1 0 0 0.1 0.28  2

#************************************************************************75
zone
1
nnum
1
1
2
nnum
1
50

#************************************************************************75
flxo
5
1 2 
3 4
20 21
48 49
49 50
#************************************************************************75
text
Time Parameters

time
1.e-2  1.e04       200    01  1995  5  0.0 
0.5 -2.0 1.0 5
1.0 -2.0 1.0 5
5.0 -2.0 1.0 5
10.0 -2.0 1.0 5

text
Numerics

ctrl
-20   1.e-4 1   40 bcgs
1   0 0 1
0 
1.0   0.0  1.
15   1.5  1.e-08 1.0e8    
0  0
iter
1.e-5 1.e-5 1.e-2 -1.e-6 1.0
00 0 0 10 14400.
sol
1    -1
#************************************************************************75
ngas
3
1 0 0 -999.


-1 0 0  0.5116e-4 0.
-2 0 0  -0.5116e-4 0.

hflx
-1 0 0 10.  1.0   
-2 0 0 25.  1.0   

adif
5.d-3
stop
</pre>

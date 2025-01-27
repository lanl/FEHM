---
title : sptr_btc
layout : page_getting-started
hero_height: is-hidden
---

# sptr_btc

 **Test sptr_btc**

This test case is constructed from the VV Test Suite sptr_btc problem which uses input files from the streamline test.
The input files are copied and used directly in this test case. These include sptr.geom, valid1.ini, and valid1.stor.

Run tests using sptr_long1alt.dat, sptr_long1btc.dat, and sptr_long1altxyz.dat.

Test Directory: [FEHM/fehmpytests/sptr_btc](https://github.com/lanl/FEHM/tree/master/fehmpytests/sptr_btc)


### Example File sptr_long1alt.dat
<pre>
*** Validation1 Test Problem: 3-D Homogeneous Flow and Transport ***
sol
    1   -1
node
5
-1 -1 -1 -1 -1
0.      0.      -10.
5000    0.      -10.
10000   0.      -10.
15000   0.      -10.
20000   0.      -10.
pres
    1   0   0   2.200    44.00   1
    0
rock
 1   0   0   2530.0  780.0e20   0.0283

perm
 1   0   0   1.0e-12  1.0e-12  1.0e-12

cond
    1   0   0  1.4e-00   1.4e-00   1.4e-00
    0
time
500000. 500000.   500   500   1997  4 0.
     0
ctrl
   10   1.e-4   12
    1    0    0    1
    0
   1.0   0    1.0
    7   2.0   0.001   1.e10
    0    +2
iter
1.e-5 1.e-5 1.e-5 1.e-9 1.1
1 0 0 0 5000.
zone
31             #Upstream Boundary
-100.      100.       100.       -100.
-100.      100.       100.       -100.
-5000.   -5000.      5000.       5000.
-5000.   -5000.      5000.       5000.
+1e12 +1e12 +1e12 +1e12
-1e12  -1e12 -1e12 -1e12
32             #Downstream Boundary
19900.      20100.       20100.       19900.
19900.      20100.       20100.       19900.
-5000.   -5000.      5000.       5000.
-5000.   -5000.      5000.       5000.
+1e12 +1e12 +1e12 +1e12
-1e12  -1e12 -1e12 -1e12

flow
     -31      0      0      2.3767     -44.00     100.00
     -32      0      0      2.0000     -44.00     100.00

zone
5
nnum
4
1938  3162  7089 51306

sptr
 1.728e8	0	0
 0.25         0
tprp
2	0.	0. 0. 0. 0. 0.

1	0	0	1

zbtc alt
1
5
 1000       1
1	100	100
0.	1.	-8.
10.	1.	-1.
1888 0.  0.1  -5.5
1888 0.  0.2  -5.5
1888 0.  0.3  -5.5
1888 0.  0.4  -5.5
1888 0.  0.5  -5.5
1888 0.  0.6  -5.5
1888 0.  0.7  -5.5
1888 0.  0.8  -5.5
1888 0.  0.9  -5.5
1888 0.  1.  -5.5
   3112         0.0      -400.1        -12.0
   3112         0.0      -400.2         -12.0
   3112         0.0      -400.3         -12.0
   3112         0.0      -400.4         -12.0
   3112         0.0      -400.5         -12.0
   3112         0.0      -400.6         -12.0
   3112         0.0      -400.7         -12.0
   3112         0.0      -400.8         -12.0
   3112         0.0      -400.9         -12.0
   3112         0.0      -401.0         -12.0
 7039         0.0       400.1       -25.0
 7039         0.0       400.2       -25.0
 7039         0.0       400.3       -25.0
 7039         0.0       400.4       -25.0
 7039         0.0       400.5       -25.0
 7039         0.0       400.6       -25.0
 7039         0.0       400.7       -25.0
 7039         0.0       400.8       -25.0
 7039         0.0       400.9       -25.0
 7039         0.0       401.0       -25.0
 51256         0.0     -2800.1      -200.0
 51256         0.0     -2800.2      -200.0
 51256         0.0     -2800.3      -200.0
 51256         0.0     -2800.4      -200.0
 51256         0.0     -2800.5      -200.0
 51256         0.0     -2800.6      -200.0
 51256         0.0     -2800.7      -200.0
 51256         0.0     -2800.8      -200.0
 51256         0.0     -2800.9      -200.0
 51256         0.0     -2801.0      -200.0

stop
</pre>

---
title : uz_test
layout : page_getting-started
hero_height: is-hidden
---

# uz_test



Test Directory: [FEHM/fehmpytests/uz_test](https://github.com/lanl/FEHM/tree/master/fehmpytests/uz_test)


### Example File uz_test.dat 
<pre>
#  france uz setup
#cont
#tec 10000 100
#1 1 0.
#head 1000.
airw
3
20.0 0.1
#........................................................
cond
  1 0 0  2. 2. 2.

rock
  1    0     0    2.50E+03     1.00E+03    0.25

#.....................................................
#wtsi
#1 20 0           #1 zone, 20 max it, 0=non-isotropic (only k(xy) is affected)
#44
#
#
zone
6
 0  1000 1000 0
 11 11 8.5 8.5
7
  0  1000 1000 0
 11 11 9.95 9.95
1
    0 .5 .5 0
   0 0 6 6

# zone 1 is the well
# zone 7 is the top face
# zone 6 is above the water table
# zone 3 is 10m from the well
# zone 4 is 30m from the well
pres
    1   0   0   .5  20.0 1
    -6   0   0   0.1  0.0 2  # uz pres macro  0.1 = initial pressure, 0.0=init sat 2=2 phase
    -7   0   0   0.1  0.0 2  # uz pres macro

flow
-7 0 0 0.1  -1.  100.  #  top face; keep at pressure of 0.1, don't allow water in

boun
file
../input/well-boun_gaz.txt
flxz
2
7 1
zonn
21
9.9 10.5 10.5 9.9
5 5 5.5 5.5

#21
# 0 150 150 0
# 5 5 10 10
#
ppor
-1
1 0 0 2.86e-3

#
#
perm
  1 0 0 -11. -11. -11.
-1 0 0 -11. -11. -9.

node
block
396 396 1
1 533 76

rlp
1 0 0 1. 1. 0. 0.

1 0 0 1

#3. 0.01 1. 5. 2.68 2 .011 van gan?
#1=linear rel perm, linear cap pressure, 0 residual liquid
# 0 residual vapor; 1 max liquid sat; 1 max vapor sat; 0. cap pressure; 0. sat. at
# which cap. pressure goes to zero.
time
0.001 10 10000 3 1999 5 0.

ctrl
 20  1.00000e-05 20 100 gmres
 1 0 0 2

 1.00000 2 1.00000
 20 1.200 1.000e-08 1.E10
 4  1
sol
 1 -1
iter
 1.00000e-05 1.00000e-05 0.100000 -1.00000e-4 1.1   # change from -8 to -6 (tolerance), .001 to .1
0       0       0       15       3600.
stop
</pre>

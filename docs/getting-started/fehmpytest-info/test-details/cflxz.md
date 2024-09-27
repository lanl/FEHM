---
title : cflxz
layout : page_getting-started
hero_height: is-hidden
---

# cflxz

**Concentration Zone Flux Test**


This test is constructed from cflxz_test in the VV Test Suite which compares FEHM results with solution from 1-D Difusion Model. FEHM writes files for Aqueous and Vapor Time vs Flux In/Out. This test case uses verified FEHM files from the VV Test which are compared using the Density variable.

Files compared include:
<pre>
cflxz_test.00001_con_node.dat  cflxz_test.00002_con_node.dat  cflxz_testAqueous_Species_001.cflx
cflxz_test.00001_sca_node.dat  cflxz_test.00002_sca_node.dat  cflxz_testVapor_Species_001.cflx
cflxz_test.00001_vec_node.dat  cflxz_test.00002_vec_node.dat
</pre>


Test Directory: [FEHM/fehmpytests/cflxz](https://github.com/lanl/FEHM/tree/master/fehmpytests/cflxz)


### Example File cflxz_test.dat 
<pre>
title: 1-d gas diffusion test
sol
1 -1
node
block
1 400 1

fdm
block
  1 400 1
  0   0 0
  1   1.0

  1    1.0
  -400 1.0

  1 1.0

zone
1
nnum
1
199
2
nnum
1
200
3
nnum
1
201

adif
0.1
air
2
20. 0.1
pres
    1   0   0  0.1 0.01   2

time
 0.001 5.  100000  50  2007 1

ctrl
30 1.e-06 30 100 gmres
1 0 0 2

1.00000 0. 1.00000
10 1.2 1e-16  .14
1  0
iter
1e-5 1e-5 1e-5 -1.00000e-5 1.2
0       0       0       5       14000.
perm
   1   0  0  1e-13 1e-13 1e-13

rock
    1   0  0  2000   1e20  0.5

hist
tecplot
mpa
deg
conc
end
cont
tec 2000 10000
t
press
conc
veloc
vapor
sat
po
perm
den
endtec
rlp
3   0.0001  1.0   3.0  3.0  2.  0.05

1 0 0 1

trac
0.0     1.0      1.e-7   1.0
0       1.0e6    1.0e6   0.0
85      1.2     1.e-3   1.0
2
 -1
0  0 0 1 1e-4  1.e-20 1.e-20 1.e-20

1 0 0 1

 1 200 1   7.85e-4

     1 200 1     0    0       -7.85e-4      0.     0.0001

 1
0  0 0 1 1e-4  1.e-20 1.e-20 1.e-20

1 0 0 1

 1 200 1   7.85e-4

     1 200 1     0    0       -7.85e-4      0.     0.0001

cflxz
2
1 3
stop

</pre>

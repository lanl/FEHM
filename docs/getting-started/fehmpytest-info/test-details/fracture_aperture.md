---
title : fracture_aperture
layout : page_getting-started
hero_height: is-hidden
---

# fracture_aperture

**Test fracture_aperture**


Comparison of 3D Wellbore Thermal Stress test results



Test Directory: [FEHM/fehmpytests/fracture_aperture](https://github.com/lanl/FEHM/tree/master/fehmpytests/fracture_aperture)


### Example File pmd_all.in 
<pre>
 ********* 3-D THERMAL STRESS_TEST  3-20-09  ******
text
3-D wellbore thermal stress vertical fracture
vertical fracture is at x = 0.0
box -200 < x < 0, 0 < y < 1000, 0 < z 500

sol
1 -1
zone
800     #z = 100-400
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
99.9 99.9 99.9 99.9 400.1 400.1 400.1 400.1
700   #fracture plane
200.1 199.999999 199.999999 200.1 200.1 199.999999 199.999999 200.1
-1.e6 -1.e6 100.1 100.1 -1.e6 -1.e6 100.1 100.1
229.9 229.9 229.9 229.9 270.1 270.1 270.1 270.1

#          read (inpt  ,   *) istrs,  istrs_coupl
#c istrs = 0 - skip stress solution
#c istrs = 1 - plain strain and 3-D (hookean) solution
#c istrs = 2 - plain stess (hookean) solution (must be 2-D)
#c idof_stress - identifies the processes (see below)
#c istrs_coupl - identifies when the stress solution is called
#c istrs_coupl = -2 beginning and end of simulation
#c istrs_coupl = -3 end of each time step
#c istrs_coupl = -1 at the end of the simulation
#c istrs_coupl = 0 just perform initial stress calculation
#c istrs_coupl = 1 at each iteration (sequential coupling)
#c istrs_coupl = 2 at each iteration (fully coupled)
#c
strs
1 -3
# sets a body force if gravity is set non zero
bodyforce
stresspor
initcalc
reldisp
permmodel
1
5 1.e-4 0. 0. 1.e4 1.e3 5.37 #
7 1.e-4 0. 0. 1.e4 1.e3 5.37 #
8 0. 0. 0. 1.e4 1.e3 5.37 0. 40. 2.

1 0 0 1
-800 0 0 4 #
-700 0 0 3 #
6001 6020 1 2

# elastic properties youngs' mod, poisson ratio
elastic
1 0 0 1.0e4 0.25 #
-700 0 0 1.0e4 0.25

# zone or zones (ja,jb,jc)  value, direction
# boundary conditions
# "distributed" means distributed to all nodes in zone
zonn
1   #x = 200
200.1 199.9999 199.9999 200.1 200.1 199.9999 199.9999 200.1
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6

stressboun
-1 0 0   0.  1

zonn
2    #x = 0
-0.1 0.1 0.1 -0.1 -0.1 0.1 0.1 -0.1
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6

stressboun #
distributed
-2 0 0 7.e6 -1

zonn
3     #y = 0
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
-0.1 -0.1 0.1 0.1 -0.1 -0.1 0.1 0.1
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6

stressboun
-3 0 0 0.  2

zonn
4     #y = 1000
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
999.9 999.9 1000.1 1000.1 999.9 999.9 1000.1 1000.1
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6

stressboun #
distributed
-4 0 0 -1.5e6 -2

zonn
5     #z = 0
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
-0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1

stressboun
-5 0 0   0.  3

zonn
6     #z = 500
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
499.9 499.9 499.9 499.9  500.1 500.1 500.1 500.1

stressboun
distributed
-6 0 0 -4.e6 -3

biot
1 0 0  1.7226e-5 1. # 1.7226e-5

# "-tolerance" is absolute residual, + is order of magnitude
#stresshis
#xyz 0. 0. 250. strx
#xyz 0. 0. 250. stry
#xyz 0. 0. 250. strz
#xyz -10. 0. 250. strx
#xyz -20. 0. 250. strx
#xyz -200. 0. 250. strx
#
tolerance
-1.e-4
stressend
cond
1 0 0 2.7e-00 2.7e-00 2.7e-00

eos     #liquid condition around 10 Mpa and 30 C
2 1 1
10 30.0 999.552 0.478 0. 0.1347 0. 0.004151154 0.2e-3 0. -1.e-15
1 30.0 0.1        0. 0.     2.0 0. 0. 1.e-5 0. 0.
#1 30.0 999.552 0.478 0. 0.1347 0. 0.004151154 0.7e-3 0. 1.e-12
ctrl
	6	1.e-04 75 100 gmres
	1	0	0	2

	1.0	3	1.0
	6	5.	1.e-10	1.e-1
	0	1
iter
1.e-6	1.e-6	0.001 -1.e-2 1.1
0 0	0	0	1.e11
time
  1e-1 0.2  2000	1	1993	02 0.0

rlp
1 0. 0. 1. 1. 0. 1.

1 0 0 1

zone
5     #z = 200-300
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
199.9 199.9 199.9 199.9 300.1 300.1 300.1 300.1

perm
1 0 0 1.d-15 1.d-15 1.d-15

rock
1 0 0 2700. 1000. 0.1

zone
10    #define zone 10 as all nodes
all
21   #x = 200
200.1 199.9999 199.9999 200.1 200.1 199.9999 199.9999 200.1
-1.e6 -1.e6 100.1 100.1 -1.e6 -1.e6 100.1 100.1
229.9 229.9 229.9 229.9 270.1 270.1 270.1 270.1

perm
-21 0 0 1.d-13 1.d-13 1.d-13

rock
-21 0 0 2700. 1000. 0.1

zone
2    #x = 0
-0.1 0.1 0.1 -0.1 -0.1 0.1 0.1 -0.1
-1.e6 -1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6
4     #y = 1000
-1.e6 1.e6 1.e6 -1.e6 -1.e6 1.e6 1.e6 -1.e6
999.9 999.9 1000.1 1000.1 999.9 999.9 1000.1 1000.1
1.e6 1.e6 1.e6 1.e6 -1.e6 -1.e6 -1.e6 -1.e6

pres
1 0 0 12. 50. 1

flow
-2 0 0 12. -50. 1.e-2
-4 0 0 12. -50. 1.e-2

zone
7    #injection point
list
200.0 0.0 250.


flow
-7 0 0   -10.  -30. 0.

weli # welltype, pwxy = xy plane, 0.05m= well radius, 8m = res.thick,1 = print
1 pwxy 0.05  4. 1

-7 0 0 1

zone
71
list
0.0 0.0 250.


node
block
-71 0 0
5874 5880 1
2000 2005 1
6000 6005 1

cont
avs 1000000 1000000
geom
liq
pres
temp
disp
stress
per
end
stop
#
</pre>

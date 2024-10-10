---
title : dispersion
layout : page_getting-started
hero_height: is-hidden
---

# dispersion
This test simulates the advection-dispersion transport of solutes in groundwater (saturated single-phase flow) with a uniform Darcy flux and linear equilibrium sorption. The time step is 0.1 day, and the total simulation time is 10 days. The 5-meter simulation domain is discretized into 500 cells. The column has an initial uniform concentration of 0 mg/l, and at the inlet, a constant concentration of 1.0 mg/l conservative tracer is continuously injected. The example input file below corresponds to Case 1, as shown in Table 2 of the write-up under the information folder. 

In addition, this test case is an example of the trac optional value NPRTTRC used to write specified time steps.
By default, all time steps are written.
Group 3 - IACCMX, DAYCM, DAYCMM, DAYCMX, NPRTTRC

Test Directory: [FEHM/fehmpytests/dispersion](https://github.com/lanl/FEHM/tree/master/fehmpytests/dispersion)


### Example File run1.dat
<pre>
title: Dipsersion test
text
See what velocity is being used to calculate dispersivity

perm
1     0 0      1.e-14  1.e-14  1.e-14

rock
   1 0  0   1400.   1000.   0.01

cond
   1 0 0   1.0  1.0  1.0

# top = 1   bottom = 2
#
pres
   1      0  0       0.2   20.    1
  0
flow
   1          2  1       0.21  -20.   1.e3
   1001    1002  1       0.20  -20.   1.e3

time
 1.e-2   10.  10000    1  1999  7    0.0
   0 0 0 0
ctrl
  -3   1.e-06  40  100 gmres
  1  0 0     2
    0    0    0    0
  1.0   0  1.
  7   1.4  1.e-9  0.01
  1    0
iter
  1.e-5 1.e-5 1.e-5 -1.e-5 1.2
  0 0 0 0  19400.
sol
 +1   -1
node
11
      2 100 200 300 400 500 600 700 800 900 1002
hist
days 50 1.0
con
sat
end
cont
avsx 50000    1.0
pres
sat
conc
vec
formatted
geom
xyz
endavsx
#-----------------------------
#-----------------------------
trac
0.0    1.5     1.e-7     1.0
0      1e10     3.65e10   0
50      1.5     1.e-6      0.001  1000
1
1
0       0.0     0.0     1.   0.0  0.1 0.1 0.1

1 0 0  1
0
1 0 0           0.0

1 2 1          -1.0   0. 1000.

#************************************************************************75
stop


</pre>

---
title : mptr
layout : page_getting-started
hero_height: is-hidden
---

# mptr

**Test Multi-Species Particle Tracking**

Compares the generated ptrk files with the old ptrk files knonw to be correct.

<pre>

Present-day mean inf ysw 8/30/2002; Update porosity and rock density gplu 5/8/20
# Particle tracking for TOUGH2 flow field
dpdp
file
../input/fehm_amr_base.dpdp  
perm
1  0  0   0.100E-14 0.100E-14 0.100E-14

rlp
1  0.  0.  1.  1.  0.  1.

1  0  0  1

rock
file
../input/fehm_amr_base.rock                             
flow

time
  365.25    7.305E6	 20000      1    1997      10   

ctrl
     -10  0.10E-03      40
       1       0       0       1
0
      1.00      3.00      1.00
       5  2.5  0.10E-09  1.82625E5 
       0       1
iter
  0.10E-04  0.10E-04  0.10E-04 -0.10E-03  0.12E+01
       0       0       0       0  0.14E+05
sol
       1      -1
rflo
air
-1
20.0  0.1
node
1
1
zone
file
../input/fehmn.zone2                            
mptr
file
../input/fehm_test_mptr.mptr
stop

</pre>
---
title : colloid_filtration
layout : page_getting-started
hero_height: is-hidden
---

# colloid_filtration

**Test Colloid Filtration**

Compares the generated ptrk files with the old ptrk files known to be correct.

Test Directory: [FEHM/fehmpytests/colloid_filtration](https://github.com/lanl/FEHM/tree/master/fehmpytests/colloid_filtration)


### Example File fehm_TSPA_base1.dat 

<pre>

Main data input file for FEHM in TSPA calculations
#
# Dual Permeability parameters
# Data source: DTN:
#
dpdp
file
../../_resources/fehm_amr_base.dpdp  
#
# perm macro is a placeholder - values not used in TSPA particle tracking runs
#
perm
1  0  0   0.100E-14 0.100E-14 0.100E-14

#
# rlp macro is a placeholder - values not used in TSPA particle tracking runs
#
rlp
1  0.  0.  1.  1.  0.  1.

1  0  0  1

#
# Rock properties - Bulk rock density and porosity
# Data source: DTN: LB0207REVUZPRP.002
# Other input in this macro (heat capacity) are placeholders - values
# not used in TSPA runs
#
rock
file
../input/fehm_amr_base.rock 
flow

#
# time macro - used to set time step and ending time of simulation
#
time
  365.25     5478750      20000      1    1997      10

#
# itfc macro used to assign values for colloid filtration at matrix interfaces
# Data source: DTN: LA0003MCG12213.002
#
itfc

1 1
13  14  -1
file
../../_resources/itfc_tsw5.txt
14  15  -1
file
../../_resources/itfc_tsw6.txt
15  16  -1
file
../../_resources/itfc_tsw7.txt
16  17  -1
file
../../_resources/itfc_tsw8.txt
17  18  -1
file
../../_resources/itfc_tsw9.txt
17  19  -1
file
../../_resources/itfc_tsw9.txt
18  20  -1
file
../../_resources/itfc_ch1.txt
19  20  -1
file
../../_resources/itfc_ch1.txt
20  26  -1
file
../../_resources/itfc_chz.txt
20  22  -1
file
../../_resources/itfc_chv.txt
26  27  -1
file
../../_resources/itfc_chz.txt
26  23  -1
file
../../_resources/itfc_chv.txt
27  28  -1
file
../../_resources/itfc_chz.txt
27  24  -1
file
../../_resources/itfc_chv.txt
28  29  -1
file
../../_resources/itfc_chz.txt
28  25  -1
file
../../_resources/itfc_chv.txt
29  30  -1
file
../../_resources/itfc_chz.txt
29  31  -1
file
../../_resources/itfc_chv.txt
18  21  -1
file
../../_resources/itfc_ch1.txt
19  21  -1
file
../../_resources/itfc_ch1.txt
21  22  -1
file
../../_resources/itfc_chv.txt
21  26  -1
file
../../_resources/itfc_chz.txt
22  23  -1
file
../../_resources/itfc_chv.txt
22  27  -1
file
../../_resources/itfc_chz.txt
23  24  -1
file
../../_resources/itfc_chv.txt
23  28  -1
file
../../_resources/itfc_chz.txt
24  25  -1
file
../../_resources/itfc_chv.txt
24  29  -1
file
../../_resources/itfc_chz.txt
25  30  -1
file
../../_resources/itfc_chz.txt
25  31  -1
file
../../_resources/itfc_chv.txt
30  32  -1
file
../../_resources/itfc_pp4.txt
31  32  -1
file
../../_resources/itfc_pp4.txt
32  33  -1
file
../../_resources/itfc_pp3.txt
33  34  -1
file
../../_resources/itfc_pp2.txt
34  35  -1
file
../../_resources/itfc_pp1.txt
35  36  -1
file
../../_resources/itfc_bf3.txt
36  37  -1
file
../../_resources/itfc_bf2.txt
16  40  -1
file
../../_resources/itfc_tsw8.txt
40  41  -1
file
../../_resources/itfc_tsw9.txt
41  42  -1
file
../../_resources/itfc_ch1.txt
42  22  -1
file
../../_resources/itfc_chv.txt
42  43  -1
file
../../_resources/itfc_chz.txt
43  27  -1
file
../../_resources/itfc_chz.txt
43  23  -1
file
../../_resources/itfc_chv.txt
28  44  -1
file
../../_resources/itfc_chz.txt
24  44  -1
file
../../_resources/itfc_chz.txt
44  45  -1
file
../../_resources/itfc_chz.txt
45  46  -1
file
../../_resources/itfc_pp4.txt
46  33  -1
file
../../_resources/itfc_pp3.txt

#
# ctrl macro used to set run control parameters
#
ctrl
     -10  0.10E-03      40
       1       0       0       1
0
      1.00      3.00      1.00
       5  1.5  0.10E-09  1.82625E6
       0       1
#
# iter macro used to set run control parameters
# these are for the FEHM flow simulations and hence are placeholders
# for the TSPA particle tracking runs
#
iter
  0.10E-04  0.10E-04  0.10E-04 -0.10E-03  0.12E+01
       0       0       0       0  0.14E+05
#
# sol macro used to set run control parameters
#
sol
       1      -1
#
# rflo macro signifies that the flow field is being read in, not conputed in FEHM
#
rflo
#
# air macro signifies that this is a two-phase problem
#
air
-1
20.0  0.1
#
# node macro used to request output at particluar nodes - node 1 is being
# requested here (placeholder)
#
node
1
1
#
# This call to the zone macro has the original layers defined (the same as
# those used in defining the layers above), followed by zones that define
# the input and output bins via the use of node lists. This overwriting
# is needed so that these zones are specified for the mptr input and
# the subsequent model run that ensues
#
zone
file
../input/fehmn.zone2_gt10
mptr
file
../input/fehm_TSPA_base1.mptr
stop

</pre>

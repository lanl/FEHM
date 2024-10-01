---
title : henrys_law
layout : page_getting-started
hero_height: is-hidden
---

# henrys_law

**Test the Henrys law Problem**

Comparison of Model and Analytical Solution for Concentration vs Time for 1) Air Movement through stagnant water 2) Water Movement through stagnant air 3) Liquid phase reaction. FEHM input files are henry1.in, henry2.in, and henry3.in

<pre>
==> henry1.in <==
***** Test Suite 1.1  ******* (old 1.2)
#      A two phase, air-water simulation was setup to verify
#    Henry's law.  A one-dimensional, plug flow reactor was
#    used for the following simulation.
#     See Bruce for full documentation
#    A Henry's law constant was chosen so that half of the species
#    resides in the vapor and half in the liquid.  In the simulation,
#    the tracer will be exchanged between the flowing vapor and the
#    stagnant liquid.  Therefore, the mean residence time of the
#    species should be approximately twice that of the vapor only
#    species run in test suite 1.1 (since all other parameters
#    were kept unchanged).  The run confirms that Henry's law is
#    correctly establishing equilibrium between the two phases.

==> henry2.in <==
***** Test Suite 2.1  ******* (old 2.2)
#     A two phase, air-water, simulation was setup to verify
#  Henry's law.  A one-dimensional, plug flow reactor was
#  used for the following simulation.
#     See Bruce for full documentation
#  A Henry's law constant was chosen so that half of the species
#  exists in the vapor and half in the liquid.  In the simulation,
#  tracer is exchanged between the flowing liquid and the
#  stagnant vapor.  Therefore, the mean residence time of the species
#  should be approximately twice that of the liquid only species run
#  in test suite 2.1 (since all other parameters are unchanged).  The
#  run confirms that Henry's law is correctly simulating equilibrium
#  between the stagnant air and the moving water.

==> henry3.in <==
***** Test Suite 3.1: TEST Henry's Law with reaction  (old 3.3)*******
#    The same two phase simulation used in test suite 2 was used
#  for the following simulation.  These runs are designed to
#  test whether the chemical reaction portion of the code works
#  for a Henry's law species.
#    See Bruce for full documentation
#  Test Suite 3.2 was run again but the reaction is run in the
#  reverse direction to confirm the reverse reaction capability
#  of the reaction module.  Test suite 3.2 and 3.3 should produce
#  identical breakthrough curves (except for round off error).
</pre>


Test Directory: [FEHM/fehmpytests/henrys_law](https://github.com/lanl/FEHM/tree/master/fehmpytests/henrys_law)


### Example File  henry1.in
<pre>
***** Test Suite 1.1  ******* (old 1.2)
#      A two phase, air-water simulation was setup to verify
#    Henry's law.  A one-dimensional, plug flow reactor was
#    used for the following simulation.
node
  1
  201
perm
  1  402 1 1.e-11 1.e-11 1.e-11
000
rock
  1  402 1    2500.   1000.   0.05
000
rlp
   1   0.3   0.3   1.0   1.0  0.0  0.0  0.0
000
  1  402 1    1
000
cond
  1  402 1  1.73   1.73   1.73
000
pres
  1  402 1   .1     0.2   2
000
airwater
   -2
   20.0   0.1
flow
  1  402 201   -1.e-4     0.   0.
201  402 201    1.e-4     0.   0.
000
sol
    1    -1
cont
  999999999  1.0
ctrl
   40   1.e-8   8
   1   402   1   1
000
   1.0   3.0   0.75
   10000   2.0   0.0000001   2000.0
1   0
time
   25.0  375.012   5000   1   9999   01
000
trac
0. 1. 1.e-6 .5
375.00001 375.012 375.00001 375.012
10000 1.2 1.e-6 1.e-6
1
-2
0 0 0 1  1.e-9 .033333 .033333 .033333
0 0 0 1  1.e-9 .033333 .033333 .033333

1 402 1 1

1 33.64 0
1 402 1 0.

1 402 201 1. 375.00001 375.01

stop
</pre>

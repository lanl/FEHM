---
title : transport3d
layout : page_getting-started
hero_height: is-hidden
---

# transport3d

**Test transport3d**


Three-Dimensional Radionuclide Transport Problem - trac_rlp.
        Comparison of FEHM and TRACRN for Concentration vs Time

Run 3 tests using 3d_trac_rlp.trc,  3d_trac_rlpm.trc, and  3d_trxn_rlp.trc for comparisons.

Test Directory: [FEHM/fehmpytests/transport3d](https://github.com/lanl/FEHM/tree/master/fehmpytests/transport3d)

Note: The VV windows and VV linux transport3d tests and perl scripts are very different. Use this test case to test both.
The VV test problem consists of 2 parts, the first generates a flow field which is used as input. Here the .ini file is read as input.

Issues: Depending on version, the .trc file formats may differ. Old versions use a single line for 6 data values, new windows uses double lines with 3 values each.


### Example File 3d_trac_rlp.dat 
<pre>
********* Tracer problem 5/16/95 gaz **********
text
3-D Tracer Problem for Verification Report

sol
	1	-1
zone
 1
	0.	100.0	100.0	0.
	0.	100.0	100.0	0.
	0.	0.	100.0	100.0
	0.	0.	100.0	100.0
	100.0	100.0	100.0	100.0
	0.	0.	0.	0.
 -1
init
	1	20.	0.	0.	0.	0.	0.	0.
air
	-2
	20.	0.1
pres
	1	0	0	1	20	1

ctrl
	-10	1.e-04	12
	1	0	0	2

	1.0	3	1.0
	9	4.0	1.e-11	1e30
	0	1
rlp
	3	0.277	1.0	3.34	1.982	2.00	0.01

	1	0	0	1

rock
	1	0	0	2700.	1e15	0.3

cond
	1	0	0	2.7e-00	2.7e-00	2.7e-00

perm
	1	0	0	1.d-12	1.d-12	1.d-12

zone
 1
	0.	100.0	100.0	0.
	0.	100.0	100.0	0.
	0.	0.	100.0	100.0
	0.	0.	100.0	100.0
	100.0	100.0	100.0	100.0
	0.	0.	0.	0.
 2
	20.0	30.0	30.0	20.0
	20.0	30.0	30.0	20.0
	20.0	20.0	30.0	30.0
	20.0	20.0	30.0	30.0
	100.0	100.0	100.0	100.0
	99.0	99.0	99.0	99.0
 3
	60.0	90.0	90.0	60.0
	60.0	90.0	90.0	60.0
	60.0	60.0	90.0	90.0
	60.0	60.0	90.0	90.0
	1.0	1.0	1.0	1.0
	0.0	0.0	0.0	0.0
 -1
flow
	-2	0	1	-1.e-4	1.00	0.0
	-3	0	1	2.0	1.0	1.

node
 6
	841	2251	4145	6039	7933	9827
time
	3.65e03	2.826250e6 	2000	1	1993	02 0.

iter
	1.e-6	1.e-6	1.0	-1.e-9	1.1
	0	0	0	0	1.e11
itup
	10
trac
	0	1	1.e-2	1.
	1e6	1e30	1.005e6	1e30
	10	2.0	1e3	3.6525e3
	3
 1
	0	0	0	1	1.e-10	5.	5.	5.

	1	0	0	1

	1	0	0	0.0

	-2	0	0	1.	1e6	1.36525e6

 1
	0	0	0	1	1.e-10	5.	5.	5.

	1	0	0	1

	1	0	0	0.0

	-2	0	0	1.	1e6	1.36535e6

 1
	0	0	0	1	1.e-10	5.	5.	5.

	1	0	0	1

	1	0	0	0.0


rxn
** NCPLX, NUMRXN
         0    1
** Coupling of the aqueous components (dRi/dUj)
3
1 0 0
0 1 0
0 0 1
** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
         1   Cons             0       0     1.e-9
         2   Am-241           0       0     1.e-9
         3   Np-237           0       0     1.e-9
** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immoblie spec, NIMM rows)
** IDVAP(IV), VAPNAM(IM), VAPPRT(IV) (ID # and name of vapor spec, NVAP rows)
** Skip nodes
0
** RSDMAX
   1.0e-10
******  Chemical reaction information  ********
** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
** CKEQ(IX) (Stability constants, NCPLX rows)
** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
** RADIOACTIVE DECAY (type 5) **
         5
** Where does the reaction take place ? ***
 1 0 0

** Half Life (years)**
         432
** Type of Reaction (solid, liquid, vapor) **
         1
** Parent and Daughter Product **
         2     3
stop

</pre>

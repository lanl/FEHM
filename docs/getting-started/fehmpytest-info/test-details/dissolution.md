---
title : dissolution
layout : page_getting-started
hero_height: is-hidden
---

# dissolution

**Test the Dissoultion Macro**

A one-dimensional transport simulation of calcite (CaC03(s)) dissolution is tested. Profiles of concentration versus reactor length, at selected times, will be compared against the analytical solution.

Details of this test are described in the FEHM V2.21 Validation Test Plan on pages 93-95 (STN: 10086-2.21-00, Rev.No. 00, Document ID: 10086-VTP-2.21-00, August 2003.md)



Test Directory: [FEHM/fehmpytests/dissolution](https://github.com/lanl/FEHM/tree/master/fehmpytests/dissolution)


### Example File dissolution.in 
<pre>
Check FEHM dissolution versus analytical solution: Matches Peter 51 node
cont
surf 10000 1.e20
geom
liquid
concentration
formatted
endavs
cond
1 102 1 2.7 2.7 2.7

ctrl
50 1e-6 8
1 102 1 2

1 0 0.5
25 1.2 1.1574e-6 10000
1 0
iter
1.0e-5 1.0e-5 1.0e-5 1.0e-9 1.1
0 0 0 0 1e9
flow
1 102 51 -0.0014992 -25. 0.
51 102 51 1. -25 -1

init
1. 25 25 0 1000 25 0 0
node
9
1 5 6 7 8 9 10 11 51
perm
1 102 1 5.0e-10 5.0e-10 5.0e-10

rock
1 102 1 1800. 1000. 0.32

sol
1 -1
time
1.e-3 4.157407407 10000 1 93 5
1.23148148 -1. 1 1
1.69444444 -1. 1 1
2.157407407 -1. 1 1

trac
0 1 1.e-7 0.5
1. 1.e20 1. 1000.
50 1.2 1.1574e-6 1.1574e-3
2
1
0 0 0 1 1.e-9 0.0067 0.0 0.0

1 0 0 1

1 0 0 6.26e-5


0
1 0 0 2e-5


rxn
** NCPLX,NUMRXN
         0, 1
** Coupling of the aqueous components (dRi/dUj)
1
1
** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
    1     Np[aq]    0      0    1.e-9
** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immobile spec, NIMM rows)
    1  Np[s]    0
** IDVAP(IV), VAPNAM(IM),VAPPRT(IV) (ID # and name of vapor species, NVAP rows)
** skip nodes? **
0
** RSDMAX
   1.0e-13
******  Chemical reaction information  ********
** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
** CKEQ(IX), HEQ(IX) (Stability constants and Enthalpys, NCPLX rows)
** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
** Precipiation/Dissolution REACTION (type 7) **
        7
** Where does the reaction take place ? ***
1 0 0

** immobile species participating in reaction **
        1
** the number of total aqueous species in reaction **
        1
** total aqueous species in reaction **
        1
** stoichiometry of the immobilie species **
        1
** stoichiometry of the aqueous species **
        1
** solubility product **
	6.26e-5
** rate constant (moles/(m^2*sec)) **
	100
** surface area of the mineral (m^2) **
	1
stop
</pre>

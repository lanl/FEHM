---
title : saltvcon
layout : page_getting-started
hero_height: is-hidden
---

# saltvcon

**Test the Salt Variable Conductivity Macro**

Tests the calculations of thermal conductivity of crushed and intact salt.

Intact salt:

``kxi = k_{t-300}(300/T)^1.14`` *Munson et al. (1990) Overtest for Simulate Defense High-Level Waste (Laboratories, SAND89-2671)*

Thermal conductivity of crushed salt from Asse mine:

``kx_asse = -270*phi^4+370*phi^3-136*phi^2+1.5*phi+5`` *Bechtold et al. (2004) Backfilling and sealing of underground respositories for radioactive waste in salt (BAMBUS II project), EUR 20621, ISBN 92-894-7767-9*

``kx = (k_{t-300}/kx_asse)*(300/T)^1.14`` if kx is less then 1.e-6, set to 1.e-6. If porosity is greater than 0.4, it is truncated as 0.4 since the Kx relationship is only valid within this range.

The excel spreadsheet /information/saltvcon.xlsx contains the associated calculations. 

Test Directory: [FEHM/fehmpytests/saltvcon](https://github.com/lanl/FEHM/tree/master/fehmpytests/saltvcon)


### Example File intact.dat 

<pre>

title: 1-d heat pipe calculation               11/28/12
#   All nodes are crushed salt
#************************************************************************75
zone
	1
nnum
	6 1 2 3 4 5 6

salt
saltvapr
1
saltnum
permavg 1.0
poravg  1.0
pormin 1.d-5

saltvcon
  3   26.85  5.4 1.14
  4   26.85  1.08 -270. 370. -136. 1.5 5  1.14

  1 0 0 1
  1 0 0 2

saltppor 
7
  1 0 0   4.866e-9 4.637 1.e-3  0.20
 
saltadif
 333
saltend    
node
6
1 2 3 4 5 6
perm 
   1 0 0 -14 -14 -14
   1 0 0 -13 -13 -13

rlp
        3  0.05 1.0  4   1.56   -10.  0.06
        1  0.0  0.   1.0    1.0   0.15  1.0

        1      0       0       2

rock  
 1 1 1   2165.  931.  0.01
 2 2 1   2165.  931.  0.1
 3 3 1   2165.  931.  0.3
 4 4 1   2165.  931.  0.5
 5 5 1   2165.  931.  0.7
 6 6 1   2165.  931.  0.9

#vcon
#  3   26.85  5.4   1.14
#  4   26.85  1.08 -270. 370. -136. 1.5 5  1.14
#
#  1 0 0 2
#
#ppor
#7
#  1 0 0   4.866e-9 4.637 1.e-3  0.20
#
flxo 
5
1 2 
2 3
3 4
4 5
5 6
# - - - - - - - - - - - - - - - - 
time
1.e-4  200.    1    01  1995  5  0.0 
0.5 -2.0 1 1
1.0 -2.0 1 1
5.0 -2.0 1 1
10.0 -2.0 1 1

ctrl
  -15   1.e-04   24  100 gmre
1   0 0 2
0 
1.0   0.0  1.
15   1.5  1.e-09 1.  
0 +1 
iter
1.e-5 1.e-5 1.e-3 -1.e-3 1.2
00 0 0 10 1000.
sol
1    -1
#- - - - - - - - - - - - - - -
#vapl    #gaz debug comment out
#adif
#333
#- - - - - - - - - - - - -
pres       #initial pres sat
1 0 0 0.1   0.10  2

ngas  reset P
3
1  1  1 -20.
2  2  1 -40.
3  3  1 -60.
4  4  1 -80.
5  5  1 -100.
6  6  1 -120.



hflx
 1 1 1 20.   1.e6
 6 6 1 120.   1.e6

#- - - - - - - - - - - - -
#cont
#avsx  100  10000.
#temp
#sat
#porosity
#perm
#mat
#conc
#pres
#vap  
#density
#liquid
#end
#cden now in moles 226/7.5 = *****************************************
cden
1
30.1
#***********************************************************
trac
 0  1  1.e-7  1.0
 0. 3652.5  1.e6  1.e6
 50  1.6 1.e-1 1.  1
 2
 1
0 0 0 1  1.e-9 .33333 .33333 .33333

  1 0 0 1

  1 0 0  6.16 


0
  1 0 0 17.2414


rxn
** NCPLX,NUMRXN
         0, 1
** Coupling of the aqueous components (dRi/dUj)
1
1
** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
    1     1    0      0    1.e-9
** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immobile spec, NIMM rows)
    1  2    0
** IDVAP(IV), VAPNAM(IM),VAPPRT(IV) (ID # and name of vapor species, NVAP rows)
** skip nodes? **
0   -1
** RSDMAX
   1.0e-9
******  Chemical reaction information  ********
** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
** CKEQ(IX), HEQ(IX) (Stability constants and Enthalpys, NCPLX rows)
** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
** Precipiation/Dissolution REACTION (type 7) **
        8
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
lookup 8
     10 40 100 150 200 250 300 350
    6.12 6.23  6.65 7.21  8.00  9.06  10.45  12.27
porosity change
** molecular weight of mineral (kg/mol), density of mineral (kg/m^3) SALT Wikipedia**
0.0558 , 2165.
** rate constant (moles/(m^2*sec)) **
        0.01
** surface area of the mineral (m^2) **
        1
stop

** rate constant (moles/(m^2*sec)) **
        0.1
        
        

</pre>

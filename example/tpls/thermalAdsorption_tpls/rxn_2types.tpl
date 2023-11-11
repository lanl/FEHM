ptf %
rxn
**[1] NCPLX NUMRXN (2 Langmuir reactions) **
0 2   
**[2] NGROUPS GROUP (Coupling of aqueous components dRi/dUj)**
1
1
**[3] IDCPNT CPNTNAM IFXCONC CPNTPRT CPNTGS **
1      methane_l    0        0 1e-9
**[4] IDCPLX CPLXNAM CPLXPRT **
**[5] IDIMM IMMNAM IMMPRT (sorbed/solid phase)** 
1      methane_s   0
**[6] IDVAP VAPNAM VAPPRT ** 
**[7] ISKIP **
0
**[8] RSDMAX **
1.e-9
******  Chemical reaction information  ********
**[9] LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
**[10] CKEQ(IX), HEQ(IX) (Stability constants and Enthalpys, NCPLX rows)
**[11] STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
**[12_1] IDRXN (Langmuir) -- ACTIVE REGOLITH **
2
**[13_1] JA JB JC Where rxn takes place (GETS OVERWRITTEN BY USERR SUBROUTINE) **
 1 0 0

**[14_1] IAQUEOUS  IMMOBILE **
1, 1
**[15_1] DISTCOEFF (custom 'userr' subroutine) **
userr
file
../macros/userr_data.dat
**[16_1] RATE (k_m)**
%rate%
**[17_1] MAXCONC **
%maxconc%
**[12_2] IDRXN (Langmuir) -- CONSTANT REGOLITH **
2
**[13_2] JA JB JC Where rxn takes place ** 
 %constant_regolith_zonen% 0 0

**[14_2] IAQUEOUS  IMMOBILE **
1, 1
**[15_2] DISTCOEFF ** 
%constant_distcoeff%
**[16_2] RATE (k_m)**
%rate%
**[17_2] MAXCONC **
%maxconc%
stop


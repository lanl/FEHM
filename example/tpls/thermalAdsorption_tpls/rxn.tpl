ptf %
rxn
**[1] NCPLX NUMRXN **
0 1   
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
**[12] IDRXN (Langmuir)**
2
**[13] JA JB JC Where rxn takes place (GETS OVERWRITTEN BY USERR SUBROUTINE) **
 1 0 0

**[14] IAQUEOUS  IMMOBILE **
1, 1
**[15] DISTCOEFF (custom 'userr' subroutine) **
userr
file
../macros/userr_data.dat
**[16] RATE (k_m)**
%rate%
**[17] MAXCONC **
%maxconc%
stop


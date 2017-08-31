ptf %
title: Radionlucide Decay and Ingrowth learning Problem
sol
1 -1
node
file
../input/nodeline.dat
airwater
2
20. 0.1
pres
  1  0  0  1.0  0.5  2

#---------------Output----------------------
hist
days 1000 0.1
tecplot
pressure
concentration
saturation
flow
end
cont
tec 1e20  10  days
liquid
pressure
velocity
concentration
endcont
#---------------Time/Ctrl-------------------
time
 0.001 20  100000000  1  2007 1  0
 0.01  -1.2 1.0 10
 0.1  -1.2 1.0 10
 1  -1.2 1.0 10
 10  -1.2 1.0 10
 20  -1.2 1.0 10
 30  -1.2 1.0 10
 40  -1.2 1.0 10
 50  -1.2 1.0 10

ctrl
30 1.e-03 30 100 gmres
1 0 0 2

1.00000 0  1.00
10 2 1e-3  1.0
1  0
iter
1e-5 1e-5 1e-3 -1e-2  2
0       0       0       5       14000.
#----------------Perm/Rock-------------------
perm
  1  0  0  1.0e-12  1.0e-12  1.0e-12

rock
  1  0  0  2000.  1010.  0.3

trac
0.0     1.0      1.e-09   1.0
0       1.0e20   1.0e20   0.0
10      1.5      1.e-5    0.01   1
3
 -2 #135iodine
0  0 0 1 1 1.00e-9  1.e-27 1.e-27 1.e-27
0  0 0 1 1 1.24e-5  1.e-27 1.e-27 1.e-27

1 0 0 1

1 1291.631 0.
1 0 0 %C0_I%


 -2 #135xenon
0  0 0 1 1 1.00e-9  1.e-27 1.e-27 1.e-27
0  0 0 1 1 1.24e-5  1.e-27 1.e-27 1.e-27

1 0 0 1

1 1291.631 0.
1 0 0 %C0_Xe%


 -2 #135cesium
0  0 0 1 1 1.00e-9  1.e-27 1.e-27 1.e-27
0  0 0 1 1 1.24e-5  1.e-27 1.e-27 1.e-27

1 0 0 1

1 1291.631 0.
1 0 0 %C0_Cs%


rxn
** NCPLX NUMRXN **
0 2
** GROUP **
1
1 1 1
** IDCPNT CPNTNAM IFXCONC CPNTPRT CPNTGS **
1 135iodine 0 0 1.e-9
2 135xenon  0 0 1.e-9
3 135cesium 0 0 1.e-9
** IDCPLX CPLXNAM CPLXPRT **
** IDIMM IMMNAM IMMPRT **
** IDVAP VAPNAM VAPPRT **
** ISKIP **
0
** RSDMAX **
1e-9 -1
** extra heading **
** Omit group 9 **
** Omit group 10 **
** Omit group 11 **
** IDRXN **
5
** JA JB JC **
 1 0 0

** HALFLIFE **
%thalf_I%
** RXNTYPE **
1
** PARENT DAUGHTER **
1 2
** IDRXN **
5
** JA JB JC **
 1 0 0

** HALFLIFE **
%thalf_Xe%
** RXNTYPE **
1
** PARENT DAUGHTER **
2 3
stop



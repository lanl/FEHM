
These files are from George running non-darcy xfehm Jun 11 2025

In summary the iterations should finish much below max 600
High Beta values should have higher max pres at end

-------------------------------
liq_ndar_10.out_no_nd ndar OFF has lowest pres:

Node   P (MPa)   P Cap (MPa) P Liq (MPa) Vapor (kg/s)     S gas        Tot (mfrac)    Liq (mfrac)
   1   3.513       0.000       3.513       0.000          0.0000        6.66837E-21     0.0000
 Net kg water discharge (total out-total in):  -0.197636E+03
 Net kg air discharge   (total out-total in):   0.000000E+00
 simulation ended: days   50.000000000000000000000     timesteps    77
 total N-R iterations =        154
 total solver iterations =        340

-------------------------------
liq_ndar_10.out_nd10 ndar 1.d10

Node   P (MPa)   P Cap (MPa) P Liq (MPa) Vapor (kg/s)     S gas        Tot (mfrac)    Liq (mfrac)
   1   16.70       0.000       16.70       0.000          0.0000        6.61884E-21     0.0000
 Net kg water discharge (total out-total in):  -0.636673E+02
 Net kg air discharge   (total out-total in):   0.000000E+00
 simulation ended: days   50.000000000000000000000     timesteps    77
 total N-R iterations =        159
 total solver iterations =        515

-------------------------------
liq_ndar_10.out_nd11 ndar 1.d11

Node   P (MPa)   P Cap (MPa) P Liq (MPa) Vapor (kg/s)     S gas        Tot (mfrac)    Liq (mfrac)
   1   127.0       0.000       127.0       0.000          0.0000        6.23155E-21     0.0000
 Net kg water discharge (total out-total in):  -0.113088E+03
 Net kg air discharge   (total out-total in):   0.791942E-16
 simulation ended: days   50.000000000000000000000     timesteps    77
 total N-R iterations =        163
 total solver iterations =        753

---
title : vapor_extraction
layout : page_getting-started
hero_height: is-hidden
---

# vapor_extraction

**Test vapor_extraction**


Comparison of Model and Analytical Solution for Vapor Pressure vs Position.
Run tests using vapextract_aniso.in and vapextract_iso.in.


**Soil Vapor Extraction -- Anisotropic Case**

Comparison of fehmn results with analytic solution obtained by
Shan, Falta and Javandel: Analytic Solutions of Steady State Gas Flow
to a Soil Vapor Extraction
Well in the Unsaturated Zone; Lawrence Berkely Lab, 06/19/91
Very fine spacing near the borehole in both r and z.
Source strength is scaled according to the variable element lengths
(see flow parameters), 0.00556 kg/s for elements of .5 m length, total flow = .05kg/s


**Soil Vapor Extraction -- Isotropic Case**

Comparison of fehmn results with analytic solution obtained by
Shan, Falta and Javandel: Analytic Solutions of Steady State Gas Flow
to a Soil Vapor Extraction
Well in the Unsaturated Zone; Lawrence Berkely Lab, 06/19/91
Very fine spacing near the borehole in both r and z
source strength is scaled according to the variable element lengths
(see flow parameters), 0.01 kg/s for elements of .5 m length.



Test Directory: [FEHM/fehmpytests/vapor_extraction](https://github.com/lanl/FEHM/tree/master/fehmpytests/vapor_extraction)


### Example File vapextract_aniso.in
<pre>
Soil Vapor Extraction -- Anisotropic Case
node
1
500
cont
avs 1000 730.
vapor
pressure
endavs
sol
   1    -1
pres
  1  1160 1  .101325  0.05  2
  1    40 1  .101325  0.05 +2

airwater
2
10.0 0.101325
rlp
3  0.10 .99 0.005 1.8 -2.00  0.1000           ! vangenutchen matrix
3  0.10 .99 0.005 1.8 2.00  0.1000           ! vangenutchen matrix

  1  1160 1  1  ! matrix

rock
1  1160  1 2700. 1000. 0.4  !  matrix

perm
  1  1160  1  1.0e-11  1.0e-12  1.0e-11
  321 801 40  1.0e-08  1.0e-08  1.0e-08

flow
  1    40 1  .101325  0.05 1.e+2
  321 321  1  0.001515  0.  0.
  361 361  1  0.001515  0.  0.
  401 401  1  0.002273  0.  0.
  441 441  1  0.004545  0.  0.
  481 481  1  0.006061  0.  0.
  521 521  1  0.006061  0.  0.
  561 561  1  0.006061  0.  0.
  601 601  1  0.006061  0.  0.
  641 641  1  0.006061  0.  0.
  681 681  1  0.004545  0.  0.
  721 721  1  0.002273  0.  0.
  761 761  1  0.001515  0.  0.
  801 801  1  0.001515  0.  0.

time
 1.0e-2   730   50    2  1992   11

ctrl
  09   1.e-05   13
  1  1160 1  2
    0    0    0    0
  1.0   2   1.0
  07   2.   1.e-11 1.0e7
  4   0
iter
1.e-5 1.e-5 1.e-5 1.e-07 1.2
0 0 0 0 40.
stop

</pre>
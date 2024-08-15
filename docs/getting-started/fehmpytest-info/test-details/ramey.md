---
title : ramey
layout : page_getting-started
hero_height: is-hidden
---

# ramey

**Test Temperature in a Wellbore Problem**

Compares the generated contour and history files to old contour and history file that are known to be correct. For the contour files, only the temperature values at time 2 are tested. For the history files, all temperature values are tested.

Test Directory: [FEHM/fehmpytests/ramey](https://github.com/lanl/FEHM/tree/master/fehmpytests/ramey)


### Example File ramey.in

<pre>

***** Temperature in a Wellbore (Ramey) *****
node
   3
   -1 -2 -3
   0. 0. 0.
   0. -1000. 0.
   0. -2000. 0.
cont
avs	10000	25.
geom
temperature
end
hist
deg
end
sol
  1 -1
init
      05.0 0. 20.0   0.03  2500.  20.0  0.03  0.
zone
   1
   0.105 40. 40. 0.105
   0. 0. -2000. -2000.
   2
   0. .1 .1 0.
   0. 0. -2000. -2000.

rock
   -1  0    0    2700.     1000.     1.e-10
   -2  0    0    2700.     1000.     1.000

cond
   -1  0    0   2.7d-00    2.7d-00   2.7d-00
   -2  0    0   2.7d-00    2.7d-00   2.7d-00

perm
   -1  0    0   1.d-20    1.d-20    1.d-20
   -2  0    0   1.d-08    1.d-08    1.d-08

time
   0.001 25. 9999  1  1990   2

ctrl
   15  1d-6  008
   -1  0  0   01
   -2  0  0   01

  1.0    0.0   1.00
   30      1.10     0.000001  1.0
    4  0
zone
   1
   0. .05 .05 0.
   0. 0. -10. -10.
   2
   0. .05 .05 0.
   -1990. -1990. -2000. -2000.

flow
   -1    0    0   -0.5    -20.00  0.0d00
   -2    0    0   00.00    -150.00  1.0d-00

stop
</pre>

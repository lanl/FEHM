---
title : Avdonin
layout : page_getting-started
hero_height: is-hidden
---

# Avdonin

**Test the Radial Heat and Mass Transfer Problem**

Compares the generated contour and history files to old contour and history files that are known to be correct. For contour files, only the temperature values at time 2 are tested. For history files, all temperature values are tested.

<pre>
	*** Avdonin fluid-to-rock heat conduction (radial flow)***
node
   1
   -1
   37.5	200.0	0.
cont
avs	1000	11574.074
geom
temperature
end
hist 
deg
end
sol
   1	-1
init
   5.	0.0	170. 0. 0. 170. 0. 0.
zone
   1 
   0.	1000.	1000.	0.   
   0.	0.	200.	200.

rock
   -1	0	0	2500.     1000.     0.20

cond
   -1	0	0	20.e-00    20.e-00   20.e-00

perm
   -1	0	0	1.e-12    1.e-12    1.e-12

time
   0.5 11574.074 2000 5 1989 10

ctrl
   25  1d-5  008
   -1	0	0	1

   1.0	0.0	1.0
   30	2.0	1.e-06	50.
   4	0
zone
   1
   0.	0.1	0.1	0.   
   0.	0.	200.	200.
   2
  999.9	1000.	1000.	999.9   
   0.	0.	200.	200.

flow
   -1	0	0	-10.000  -160.0  0.00
   -2	0	0	5.0      -170.0  1.e2

stop
</pre>
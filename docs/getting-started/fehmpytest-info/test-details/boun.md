---
title : boun
layout : page_getting-started
hero_height: is-hidden
---

# boun

**Test the Boundry Functionality**

Compares the generated contour files to old contour files that are known to be correct. Only the pressure and hydraulic head values at time 2 are tested.

Test Directory: [FEHM/fehmpytests/boun](https://github.com/lanl/FEHM/tree/master/fehmpytests/boun)


### Example File prob_well_boun1.dat 

<pre>

***Compare wellbore module***
cont
avs 1.e6  50.
xyz
geom
pres
head
liquid
formatted
endavs
head 0.
zone
50
all

#wtsi -1
#1 50 0.01 0.00000 1. -999. 0.0001 
air
-2
20. 0.1
pres 
     1     0   0   250.   20.   1 
         
time
1.e-3 1000 20  1 2007 10000 0.00000

iter
 1.e-06 1.e-06 1.e-2 -1.e-3 1.15
13       0       0       5       1200.   
perm
        1  0 0     -13 -13 -13  
        
rock
      1  0  0 2650. 750. 0.15

cond
  1   0     0     1 1 1     

zonn
11
-0.1 0.1 0.1 -0.1 -0.1 0.1 0.1 -0.1
-0.1 -0.1 1000.1 1000.1 -0.1 -0.1 1000.1 1000.1
-0.1 -0.1 -0.1 -0.1 30.1 30.1 30.1 30.1
12
999.9 1000.1 1000.1 999.9 999.9 1000.1 1000.1 999.9
-0.1 -0.1 1000.1 1000.1 -0.1 -0.1 1000.1 1000.1
-0.1 -0.1 -0.1 -0.1 30.1 30.1 30.1 30.1
13
-0.1 1000.1 1000.1 -0.1 -0.1 1000.1 1000.1 -0.1
-0.1 -0.1 0.1 0.1 -0.1 -0.1 0.1 0.1
-0.1 -0.1 -0.1 -0.1 30.1 30.1 30.1 30.1
14
-0.1 1000.1 1000.1 -0.1 -0.1 1000.1 1000.1 -0.1
999.9  999.9 1000.1 1000.1 999.9 999.9 1000.1 1000.1
-0.1 -0.1 -0.1 -0.1 30.1 30.1 30.1 30.1
21
list
300. 500. 100.


boun
model
ti
2 0. 1.e20
dsw
-10. -10.

-21 0 0 1

flow
 -12   0 0  200 1. 1.e-4
 -13   0 0  200 1. 1.e-4
 -14   0 0  200 1. 1.e-4
 -11   0 0  200 1. 1.e-4

#boun
#model
#ti
#2 0. 1.e20
#hd
#200. 200.
#
#-11 0 0 1
#
node
block
-21 0 0 
-11 0 0 

hist
years 1 1.e20
mpa
deg
end
ctrl
 -20 1.e-06 40 100 gmre
  1 0 0 2

 1.0 3 1.0
 20 1.500 1.0e-7 100. 
 0  0 
stop
fdm
poin
21 21 16
1 0.0
-21 1000.

1 0.0
-21 1000.

  1  200.
 -16   0.
</pre>

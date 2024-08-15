---
title : richards
layout : page_getting-started
hero_height: is-hidden
---

# richards

**Richards equation test**

Comparison of Richards equation with 2-phase solution Node by node comparison.

Test Directory: [FEHM/fehmpytests/richards](https://github.com/lanl/FEHM/tree/master/fehmpytests/richards)


### Example File rich_rlp.dat 

<pre>

NER_DATA March 2008 
# 11 unsaturated units
#01_Qal
#02_Qbt_1g Tshirege
#03_Qct
#04_Qbo
#05_Qbog
#06_Tpf (upper) 
#07_Tb4 tholeitic
#08_Tb4 alkalic
#09_Tb4 tephra
#10_Ta
#11_Tpf (lower)
#
#Air water statement
airwater
3 
20.0 0.10
#end air water statement
# 100 meter by 322m vadose zone (2-D)
zone
1  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
0. 0. 0. 0. -10 -10 -10 -10
2  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
 -10 -10 -10 -10  -45. -45. -45. -45.
3  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
 -45. -45. -45. -45. -55. -55. -55. -55.
4  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
 -55. -55. -55. -55. -115 -115. -115. -115.
5  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-105 -105. -105. -105. -125. -125. -125. -125.
6  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-115. -115. -115. -115. -145. -145. -145. -145. 
7  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-135. -135. -135. -135. -195. -195. -195. -195. 
8  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-185. -185. -185. -185. -225. -225. -225. -225. 
9  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-215. -215. -215. -215. -227. -227. -227. -227.
10  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-217. -217. -217. -217. -237. -237. -237. -237. 
11  
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-227. -227. -227. -227. -322. -322. -322. -322. 

#begin finite difference grid statement
fdm
block
25 1 322
0. 0. 0.
1 4.
-25 4.
   
1 1.

1 1.
-322 1.

#end finite difference grid statement
#
#
#
#
#begin intrinsic permeability statement
perm
1 0 0 -12 -12 -12
-1  0 0 4.62E-13      4.62E-13    4.62E-13   #kay ref Apendix A
-2	0	0	3.68e-13	3.68e-13	3.68e-13
-3  	0	0	8.82e-13	8.82e-13	8.82e-13
-4  	0	0	7.25e-13	7.25e-13	7.25e-13
-5	0	0	1.53e-13	1.53e-13	1.53e-13
-6  	0	0	4.73e-12	4.73e-12	4.73e-12
-7  	0	0	2.96e-13	2.96e-13	2.96e-13  #7,8,9 one unit in LANL model
-8 	0	0	2.96e-13	2.96e-13	2.96e-13
-9  	0	0	2.96e-13	2.96e-13	2.96e-13
-10     0       0       1.1489E-13	1.1489E-13      4.4865E-15   # leave as NER values
-11 	0	0	4.73e-12	4.73e-12	4.73e-12

#end intrinsic permeability statement
#
#specify rock or sand properties
rock  # porosoty models, use  as NER values
1 0 0 1600 1000 0.2
-1 0 0 1600 1000 0.43
-2 0 0 1600 1000 0.48
-3 0 0 1600 1000 0.53   #placeholder (where did porosity come from?)
-4 0 0 1600 1000 0.455
-5 0 0 1600 1000 0.5
-6 0 0 1600 1000 0.207
-7 0 0 1600 1000 0.1
-8 0 0 1600 1000 0.086
-9 0 0 1600 1000 0.274
-10 0 0 1600 1000 0.227
-11 0 0 1600 1000 0.1

#end rock properties
#
#
#start relative permeability statment and capillary pressure
rlp
3       0.087   1.      3.85    1.558   5       0.089  Estimate from Kay
3	0.018	1.	2.22	1.592	2.	0.019	Material 14
3	0.01	1.	1.52	1.506	2.	0.011	Material 12
3	0.026	1.	0.66	1.711	2.	0.027	Material 11
3	0.01	1.	0.081	4.026	2.	0.011	Material 10
3	0.01	1.	5.	2.68	2.	0.011	Material 6
3	0.03	1.	5.	1.5	2.	0.040	Material 7 
3	0.03	1.	5.	1.5	2.	0.040	Material 7 
3	0.03	1.	5.	1.5	2.	0.040	Material 7
3 0.001   1. 0.248418  1.689 5 0.002    #left 10 with NER data
3	0.01	1.	5.	2.68	2.	0.011	Material 6

 -1  0 0  1 
 -2  0 0  2 
 -3  0 0  3 
 -4  0 0  4 
 -5  0 0  5 
 -6  0 0  6 
 -7  0 0  7 
 -8  0 0  8 
 -9  0 0  9 
 -10 0 0  10
 -11 0 0  11

#end of intrinsic permeability statment
#
#input inital conditions
#end input initial conditions
#
zone   #define top boundary (full boundary)
21
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-1. -1. -1. -1. 0. 0. 0. 0.
22    #define top boundary (subset for recharge)
46. 54. 54. 46. 46. 54. 54. 46.
0. 0. 1. 1. 0. 0. 1. 1.
-1. -1. -1. -1. 0. 0. 0. 0.
31    #define bottom boundary (for disharge)
0. 100. 100. 0. 0. 100. 100. 0. 
0. 0. 1. 1. 0. 0. 1. 1.
-321. -321. -321. -321. -322. -322. -322. -322. 

#
#input boundary conditions 
#modify for positive head
boun
model 1
ti
2 0. 1.e20
dsw
-1.60e-3 -1.60e-3

-22 0 0 1

pres
1 0 0 0.1 .2 2
-31 0 0 0.1 0.9 2

flow #air pressure = 0.1 on top boundary
-21 0 0 0.1 -0.1 1.e-12

flow
-31 0 0 0.1  1 1.e0 

#end boundary conditions
#
#solver control options#
ctrl
10 1E-10 20 150 gmres
1 0 0 2

1 3 1 
10 1.5 1e-8 3650.E6
0 0
#
#end solver control options
#
#
#begin iter statement
iter
1.e-10 1.e-10 0.1 -1.e-3 1.2
0 0 0 20 1000
#end iter statement
#
#
#time statement
time
1e-1 9131.25 100 100 08 01 0

#end time statement
#
#node statement
#see page 8100 for explanation
node
bloc
13 8050 25 

#end node statement
hist
tecplot
mpa
saturation
end
#begin contour statement
cont
avs 10000  1.e20
l 
p
s
endavs
#end contour statement
stop
zone
50
-1.e5 1.e5 1.e5 -1.e5 -1.e5 1.e5 1.e5 -1.e5
-1.e5 -1.e5 1.e5 1.e5 -1.e5 -1.e5 1.e5 1.e5
-1.e5 -1.e5 -1.e5 -1.e5
1.e5 1.e5 1.e5 1.e5 

zone
22    #define top boundary (subset for recharge)
46. 54. 54. 46. 46. 54. 54. 46.
0. 0. 1. 1. 0. 0. 1. 1.
-1. -1. -1. -1. 0. 0. 0. 0.

trac
0 1 1e-5 1
1. 1.e8 1.e8 1.e8
50 1.2 0.1 500.
1
1
1 0 0 1 1e-15 1.0 1.0 1.0

1 0 0 1

1 0 0 0

-22 0 0  -1.0 1. 1.e20

#end tracer statement
stop
</pre>

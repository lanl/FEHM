---
title : dryout
layout : page_getting-started
hero_height: is-hidden
---

# dryout

**Test Dry-Out of a Partially Saturated Medium**

Compares the generated contour files to old contour files known to be correct. The saturation is tested for all times.

Test Directory: [FEHM/fehmpytests/dryout](https://github.com/lanl/FEHM/tree/master/fehmpytests/dryout)


### Example File dryout1.in

<pre>

Dry out medium using ngas to flow air through stagnant water
cond
1 0 0 1.7 1.7 1.7

cont
avs 1000000 1000000
geom
liquid
saturation
endavs
ctrl
10 1e-6 8
1 0 0 1

1 0 1.0
25 1.2 1.5e-8 1.5
1 0
node
5
1 51 101 151 201
flow
1 202 201 0. 1.e-20 0.

perm
1 0 0 1.e-14 1.e-14 1.e-14

rlp
1 0.3 0.3 1. 1. 93.6 100.

1 0 0 1

rock
1 0 0 2500 1000 0.05

pres
1 0 0 .1 .2 2
201 402 201 .1 .2 -2

hflx
1 202 201 20. 1.
201 402 201 20. 1.

ngas
3
1 0 0 -20.


1 202 201 -1.e-6  0.

sol
1 -1
iter
1.e-6 1.e-6 1.e-4 -1.e-3 1.1
0 0 0 0 10000000
itup
2
time
1.e-2 500. 1000 1 95 03
100. 1.5 1 1
200. 1.5 1 1
300. 1.5 1 1
400. 1.5 1 1

stop

</pre>

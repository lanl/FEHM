---
title : cellbased
layout : page_getting-started
hero_height: is-hidden
---

# cellbased

**Test the Cell-Based Particle Tracking Model**

Compares the generated output files with the old files known be correct. All values are tested for a root mean square difference of less than 0.05.

Test Directory: [FEHM/fehmpytests/cellbased](https://github.com/lanl/FEHM/tree/master/fehmpytests/cellbased)


### Example File chain.in

<pre>

Problem 1: Check of FEHMN dispersion particle tracking - Pe = 20
rest
old

cond
1 0 0 2.7 2.7 2.7

ctrl
50 1e-6 8
1 0 0 2

1 0 0.5
25 1.2 1.5e-8 1.e20
1 0
flow
1 13 12 -11.218809 -25 0
12 24 12 1. -25 -1

init
1. 25 25 0 1000 25 0 0
node
11
1 3 4 5 6 7 8 9 10 11 12
perm
1 0 0 5.0e-11 5.0e-30 5.0e-30

rock
1 0 0 2500 1000 0.3

sol
1 -1
time
2.e-3 2.e-3 1000 100 95 03 0.

ptrk
100000 123456
0. 2.e-3 1.e30 1.e30
1 0 2 2
2 0. 0.01 0.01 0.01 1.e-30 1. 0. 0.

1 0 0 1

1 1 1  1. 0. 1.e-8

stop

</pre>

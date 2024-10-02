---
title : rad_decay
layout : page_getting-started
hero_height: is-hidden
---

# rad_decay

**Test radioactive decay option in rxn macro**


The simulation is a batch reactor without flow and comparison is made to the Bateman equation. 
Decay of 135I->135Xe->135Cs is modeled.
The test ensures that FEHM matches the bateman equation within 10% relative error for all concentrations greater than 1e-6 moles/kg vapor.

H. Bateman. "Solution of a System of Differential Equations Occurring in the
Theory of Radio-active Transformations," Proc. Cambridge Phil. Soc. IS,
423 (1910) https://archive.org/details/cbarchive_122715_solutionofasystemofdifferentia1843

This test case is not party of the VV Test Suite and does not use files in /compare as reference, instead the results are compared against equation values in bateman.py


Test Directory: [FEHM/fehmpytests/rad_decay](https://github.com/lanl/FEHM/tree/master/fehmpytests/rad_decay)



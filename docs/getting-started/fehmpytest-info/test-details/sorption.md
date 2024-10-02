---
title : sorption
layout : page_getting-started
hero_height: is-hidden
---

# sorption

**Test One Dimensional Reactive Solute Transport**


This test case is constructed from the VV Test Suite sorption problem. Compares the generated tracer files to old tracer files known to be correct. All concentraction values are tested.
There are 2 runs using sorption_trac.in and sorption_trxn.in

One Dimensional Reactive Solute Transport - trac macro Test Cases:

Isotherm: Conservative     
Isotherm: Linear         
Isotherm: Langmuir        
Isotherm: Freundlich            
Isotherm: Modified Freundlich     

One Dimensional Reactive Solute Transport - trxn macro Test Cases:

Isotherm: Conservative   
Isotherm: Linear    
Isotherm: Langmuir  
Isotherm: Freundlich     
Isotherm: Modified Freundlich      

Files to compare:
<pre>
sorption_trac_Aqueous_Species_001.trc  sorption_trxn_Conservative.trc
sorption_trac_Aqueous_Species_002.trc  sorption_trxn_Freundlich.trc
sorption_trac_Aqueous_Species_003.trc  sorption_trxn_Langmuir.trc
sorption_trac_Aqueous_Species_004.trc  sorption_trxn_Linear.trc
sorption_trac_Aqueous_Species_005.trc  sorption_trxn_Mod-Freundlich.trc  
</pre>

Test Directory: [FEHM/fehmpytests/sorption](https://github.com/lanl/FEHM/tree/master/fehmpytests/sorption)



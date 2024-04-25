---
title : Macros
layout : page_getting-started
permalink: /program-specification/macros
---


# Macros by Problem Type


## Heat Conduction

| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| [title](Macros/MacroTitle.html)  | [cont](Macros/MacroCont.html)                                    |
| [boun](Macros/MacroBoun.html) or [flow](Macros/MacroFlow.html) or [hflx](Macros/MacroHflx.html)   | [finv](Macros/MacroFinv.html)                                    |
| [cond](Macros/MacroCond.html)                                                                         | [flo2](Macros/MacroFlo2.html)                                    |
| [coor](Macros/MacroCoor.html)                                                                         | [flxo](Macros/MacroFlxo.html) or [flxz](Macros/MacroFlxz.html) |
| [ctrl](Macros/MacroCtrl.html)                                                                         | [iter](Macros/MacroIter.html)                                    |
| [elem](Macros/MacroElem.html)                                                                         | [node](Macros/MacroNode.html) or [nod2](Macros/MacroNod2.html) |
| [init](Macros/MacroInit.html) or [pres](Macros/MacroPres.html)                                      | [renu](Macros/MacroRenu.html)                                    |
| [rock](Macros/MacroRock.html)                                                                         | [rflx](Macros/MacroRflx.html)                                    |
| [sol](Macros/MacroSol.html)                                                                           | text or comments (#)                                               |
| [time](Macros/MacroTime.html)                                                                         | [user](Macros/MacroUser.html)                                    |
| [stop](Macros/MacroStop.html)                                                                         | [vcon](Macros/MacroVcon.html)                                    |
|                                                                                                         | [zone](Macros/MacroZone.html) or [zonn](Macros/MacroZonn.html) |


## Water / Water Vapor / Heat Equivalent Continuum, Dual Porosity,Dual Permeability


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| [title](Macros/MacroTitle.html)                                                                              | [cden](Macros/MacroCden.html)                                       |
| [boun](Macros/MacroBoun.html) or [flow](Macros/MacroFlow.html) or [hflx](Macros/MacroHflx.html)          | [cont](Macros/MacroCont.html)                                       |
| [cond](Macros/MacroCond.html)                                                                                | [eos](Macros/MacroEos.html)                                         |
| [coor](Macros/MacroCoor.html)                                                                                | [exrl](Macros/MacroExrl.html)                                       |
| [ctrl](Macros/MacroCtrl.html)                                                                                | [finv](Macros/MacroFinv.html)                                       |
| [elem](Macros/MacroElem.html)                                                                                | [flo2](Macros/MacroFlo2.html)                                       |
| [init](Macros/MacroInit.html) or [pres](Macros/MacroPres.html)                                             | [flxo](Macros/MacroFlxo.html) or [flxz](Macros/MacroFlxz.html)    |
| [perm](Macros/MacroPerm.html)                                                                                | [fper](Macros/MacroFper.html)                                       |
| [rlp](Macros/MacroRlp.html)                                                                                  | [gdpm](Macros/MacroGdpm.html)                                       |
| [rock](Macros/MacroRock.html)                                                                                | [hflx](Macros/MacroHflx.html)                                       |
| [sol](Macros/MacroSol.html)                                                                                  | [iter](Macros/MacroIter.html)                                       |
| [time](Macros/MacroTime.html)                                                                                | [node](Macros/MacroNode.html) or [nod2](Macros/MacroNod2.html)    |
| [stop](Macros/MacroStop.html)                                                                                | [ppor](Macros/MacroPpor.html)                                       |
|                                                                                                                | [renu](Macros/MacroRenu.html)                                       |
| dual (* only)                                                                                                  | [rflx](Macros/MacroRflx.html)                                       |
| dpdp (** only)                                                                                                 | [rxn](Macros/MacroRxn.html)                                         |
|                                                                                                                | text or comments (#)                                                  |
|                                                                                                                | [trac](Macros/MacroTrac.html)                                       |
|                                                                                                                | [user](Macros/MacroUser.html) or [userc](Macros/MacroUserc.html)  |
|                                                                                                                | [vcon](Macros/MacroVcon.html)                                       |
|                                                                                                                | [velo](Macros/MacroVelo.html)                                       |
|                                                                                                                | [zone](Macros/MacroZone.html) or [zonn](Macros/MacroZonn.html)    |


## Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability 


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| [title](Macros/MacroTitle.html)                                                                     | [adif](Macros/MacroAdif.html)                                       |
| [boun](Macros/MacroBoun.html) or [flow](Macros/MacroFlow.html) or [hflx](Macros/MacroHflx.html) | [cden](Macros/MacroCden.html)                                       |
| [cond](Macros/MacroCond.html)                                                                       | [cont](Macros/MacroCont.html)                                       |
| [coor](Macros/MacroCoor.html)                                                                       | [eos](Macros/MacroEos.html)                                         |
| [ctrl](Macros/MacroCtrl.html)                                                                       | [finv](Macros/MacroFinv.html)                                       |
| [elem](Macros/MacroElem.html)                                                                       | [flo2](Macros/MacroFlo2.html)                                       |
| [init](Macros/MacroInit.html) or [pres](Macros/MacroPres.html)                                    | [flxo](Macros/MacroFlxo.html)                                       |
| [ngas](Macros/MacroNgas.html)                                                                       | [fper](Macros/MacroFper.html)                                       |
| [perm](Macros/MacroPerm.html)                                                                       | [gdpm](Macros/MacroGdpm.html)                                       |
| [rlp](Macros/MacroRlp.html)                                                                         | [iter](Macros/MacroIter.html)                                       |
| [rock](Macros/MacroRock.html)                                                                       | [node](Macros/MacroNode.html) or [nod2](Macros/MacroNod2.html)    |
| [sol](Macros/MacroSol.html)                                                                         | [ppor](Macros/MacroPpor.html)                                       |
| [time](Macros/MacroTime.html)                                                                       | [renu](Macros/MacroRenu.html)                                       |
| [stop](Macros/MacroStop.html)                                                                       | [rflx](Macros/MacroRflx.html)                                       |
|                                                                                                       | [rxn](Macros/MacroRxn.html)                                         |
| dual (*only)                                                                                          | [szna](Macros/MacroSzna.html)                                       |
| dpdp (**only)                                                                                         | text or comments (#)                                                  |
|                                                                                                       | [trac](Macros/MacroTrac.html)                                       |
|                                                                                                       | [user](Macros/MacroUser.html) or [userc](Macros/MacroUserc.html)  |
|                                                                                                       | [vapl](Macros/MacroVapl.html)                                       |
|                                                                                                       | [vcon](Macros/MacroVcon.html)                                       |
|                                                                                                       | [velo](Macros/MacroVelo.html)                                       |
|                                                                                                       | [zone](Macros/MacroZone.html) or [zonn](Macros/MacroZone.html)    |



## Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| [title](Macros/MacroTitle.html)                                    | [bous](Macros/MacroBous.html)                                       |
| [airwater](Macros/MacroAirwater.html)                              | [cont](Macros/MacroCont.html)                                       |
| [boun](Macros/MacroBoun.html) or [flow](Macros/MacroFlow.html)   | [eos](Macros/MacroEos.html)                                         |
| [coor](Macros/MacroCoor.html)                                      | [exri](Macros/MacroExri.html)                                       |
| [ctrl](Macros/MacroCtrl.html)                                      | [finv](Macros/MacroFinv.html)                                       |
| [elem](Macros/MacroElem.html)                                      | [flo2](Macros/MacroFlo2.html)                                       |
| [init](Macros/MacroInit.html) or [pres](Macros/MacroPres.html)   | [flxo](Macros/MacroFlxo.html)                                       |
| [node](Macros/MacroNode.html) or [nod2](Macros/MacroNod2.html)   | [fper](Macros/MacroFper.html)                                       |
| [perm](Macros/MacroPerm.html)                                      | [gdpm](Macros/MacroGdpm.html)                                       |
| [rock](Macros/MacroRock.html)                                      | [head](Macros/MacroHead.html)                                       |
| [sol](Macros/MacroSol.html)                                        | [iter](Macros/MacroIter.html)                                       |
| [time](Macros/MacroTime.html)                                      | [ppor](Macros/MacroPpor.html)                                       |
| [stop](Macros/MacroStop.html)                                      | [pres](Macros/MacroPres.html)                                       |
|                                                                      | [renu](Macros/MacroRenu.html)                                       |
| dual (*only)                                                         | [rlp](Macros/MacroRlp.html)                                         |
| dpdp (only)                                                          | [rxn](Macros/MacroRxn.html)                                         |
|                                                                      | text or comments                                                      |
|                                                                      | [trac](Macros/MacroTrac.html)                                       |
|                                                                      | [user](Macros/MacroUser.html) or [userc](Macros/MacroUserc.html)  |
|                                                                      | [vapl](Macros/MacroVapl.html)                                       |
|                                                                      | [velo](Macros/MacroVelo.html)                                       |
|                                                                      | [zone](Macros/MacroZone.html) or [zonn](Macros/MacroZonn.html)    |



### Required Macros


   * [anpe](Macros/MacroAnpe.rst)
   * [boun](Macros/MacroBoun.rst)
   * [cond](Macros/MacroCond.rst)
   * [coor](Macros/MacroCoor.rst)
   * [elem](Macros/MacroElem.rst)
   * [fdm](Macros/MacroFdm.rst)
   * [flow](Macros/MacroFlow.rst)
   * [ftsc](Macros/MacroFtsc.rst)
   * [hcon](Macros/MacroHcon.rst)
   * [hyco](Macros/MacroHyco.rst)
   * [init](Macros/MacroInit.rst)
   * [ivfc](Macros/MacroIvfc.rst)
   * [perm](Macros/MacroPerm.rst)
   * [pres](Macros/MacroPres.rst)
   * [rock](Macros/MacroRock.rst)
   * [stop](Macros/MacroStop.rst)
   * [time](Macros/MacroTime.rst)



### Optional Macros


   * [adif](Macros/MacroAdif.rst)
   * [airwater or air](Macros/MacroAirwater.rst)
   * [alta](Macros/MacroAlti.rst) (Deprecated)
   * [bous](Macros/MacroBous.rst)
   * [carb](Macros/MacroCarb.rst)
   * [cden](Macros/MacroCden.rst)
   * [cflx](Macros/MacroCflx.rst)
   * [cgdp](Macros/MacroCgdp.rst)
   * [chea](Macros/MacroChea.rst)
   * [conn](Macros/MacroConn.rst)
   * [cont](Macros/MacroCont.rst)
   * [conv](Macros/MacroConv.rst)
   * [ctrl](Macros/MacroCtrl.rst)
   * [dpdp](Macros/MacroDpdp.rst)
   * [dual](Macros/MacroDual.rst)
   * [dvel](Macros/MacroDvel.rst)
   * [eos](Macros/MacroEos.rst)
   * [evap](Macros/MacroEvap.rst)
   * [exrl](Macros/MacroExrl.rst)
   * [finv](Macros/MacroFinv.rst)
   * [flgh](Macros/MacroFlgh.rst)
   * [flo2](Macros/MacroFlo2.rst)
   * [flo3](Macros/MacroFlo3.rst)
   * [floa](Macros/MacroFloa.rst)
   * [flwt](Macros/MacroFlwt.rst)
   * [flxn](Macros/MacroFlxn.rst)
   * [flxo](Macros/MacroFlxo.rst)
   * [flxz](Macros/MacroFlxz.rst)
   * [fper](Macros/MacroFper.rst)
   * [frlp](Macros/MacroFrlp.rst)
   * [gdkm](Macros/MacroGdkm.rst)
   * [gdpm](Macros/MacroGdpm.rst)
   * [grad](Macros/MacroGrad.rst)
   * [head](Macros/MacroHead.rst)
   * [hflx](Macros/MacroHflx.rst)
   * [hist](Macros/MacroHist.rst)
   * [ice](Macros/MacroIce.rst) or [meth](Macros/MacroIce.rst)
   * [intg](Macros/MacroIntg.rst)
   * [imex](Macros/MacroImex.rst)
   * [impf](Macros/MacroImpf.rst)
   * [isot](Macros/MacroIsot.rst)
   * [iter](Macros/MacroIter.rst)
   * [itfc](Macros/MacroItfc.rst)
   * [ittm](Macros/MacroIttm.rst)
   * [itup](Macros/MacroItup.rst)
   * [iupk](Macros/MacroIupk.rst)
   * [mdnode](Macros/MacroMdnode.rst)
   * [mptr](Macros/MacroMptr.rst)
   * [nfinv](Macros/MacroNfinv.rst)
   * [ngas](Macros/MacroNgas.rst)
   * [nobr](Macros/MacroNobr.rst)
   * [nod2](Macros/MacroNod2.rst)
   * [nod3](Macros/MacroNod3.rst)
   * [node](Macros/MacroNode.rst)
   * [nrst](Macros/MacroNrst.rst)
   * [para](Macros/MacroPara.rst)
   * [pest](Macros/MacroPest.rst)
   * [phys](Macros/MacroPhys.rst)
   * [ppor](Macros/MacroPpor.rst)
   * [ptrk](Macros/MacroPtrk.rst)
   * [renu](Macros/MacroRenu.rst)
   * [rest](Macros/MacroRest.rst)
   * [rflo](Macros/MacroRflo.rst)
   * [rflx](Macros/MacroRflx.rst) (Deprecated)
   * [rich](Macros/MacroRich.rst)
   * [rive](Macros/MacroRive.rst) or [well](Macros/MacroRive.rst)
   * [rlp](Macros/MacroRlp.rst)
   * [rlpm](Macros/MacroRlpm.rst)
   * [rxn](Macros/MacroRxn.rst)
   * [sol](Macros/MacroSol.rst)
   * [sptr](Macros/MacroSptr.rst)
   * [stea](Macros/MacroStea.rst)
   * [strs](Macros/MacroStrs.rst)
   * [subm](Macros/MacroSubm.rst)
   * [svar](Macros/MacroSvar.rst)
   * [szna](Macros/MacroSzna.rst)
   * [text](Macros/MacroText.rst)
   * [thic](Macros/MacroThic.rst)
   * [trac](Macros/MacroTrac.rst)
   * [trxn](Macros/MacroTrxn.rst)
   * [user](Macros/MacroUser.rst)
   * [vap](Macros/MacroVapl.rst)
   * [abou](Macros/MacroVbou.rst)
   * [vcon](Macros/MacroVcon.rst)
   * [weli](Macros/MacroWeli.rst)
   * [wflo](Macros/MacroWflo.rst)
   * [wgtu](Macros/MacroWgtu.rst)
   * [wtsi](Macros/MacroWtsi.rst)
   * [zeol](Macros/MacroZeol.rst)
   * [zneg](Macros/MacroZneg.rst)
   * [zone](Macros/MacroZone.rst)
   * [zonn](Macros/MacroZonn.rst)
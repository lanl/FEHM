---
title : Macros
layout : page_getting-started
permalink: /program-specification/macros
hero_height: is-hidden
---


# Macros by Problem Type


## Heat Conduction

| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| title  | [cont](../Macros/MacroCont.md)                                    |
| [boun](../Macros/MacroBoun.md) or [flow](../Macros/MacroFlow.md) or [hflx](../Macros/MacroHflx.md)   | [finv](../Macros/MacroFinv.md)                                    |
| [cond](../Macros/MacroCond.md)                                                                         | [flo2](../Macros/MacroFlo2.md)                                    |
| [coor](../Macros/MacroCoor.md)                                                                         | [flxo](../Macros/MacroFlxo.md) or [flxz](../Macros/MacroFlxz.md) |
| [ctrl](../Macros/MacroCtrl.md)                                                                         | [iter](../Macros/MacroIter.md)                                    |
| [elem](../Macros/MacroElem.md)                                                                         | [node](../Macros/MacroNode.md) or [nod2](../Macros/MacroNod2.md) |
| [init](../Macros/MacroInit.md) or [pres](../Macros/MacroPres.md)                                      | [renu](../Macros/MacroRenu.md)                                    |
| [rock](../Macros/MacroRock.md)                                                                         | [rflx](../Macros/MacroRflx.md)                                    |
| [sol](../Macros/MacroSol.md)                                                                           | text or comments (#)                                               |
| [time](../Macros/MacroTime.md)                                                                         | [user](../Macros/MacroUser.md)                                    |
| [stop](../Macros/MacroStop.md)                                                                         | [vcon](../Macros/MacroVcon.md)                                    |
|                                                                                                         | [zone](../Macros/MacroZone.md) or [zonn](../Macros/MacroZonn.md) |


## Water / Water Vapor / Heat Equivalent Continuum, Dual Porosity,Dual Permeability


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| title                                                                              | [cden](../Macros/MacroCden.md)                                       |
| [boun](../Macros/MacroBoun.md) or [flow](../Macros/MacroFlow.md) or [hflx](../Macros/MacroHflx.md)          | [cont](../Macros/MacroCont.md)                                       |
| [cond](../Macros/MacroCond.md)                                                                                | [eos](../Macros/MacroEos.md)                                         |
| [coor](../Macros/MacroCoor.md)                                                                                | [exrl](../Macros/MacroExrl.md)                                       |
| [ctrl](../Macros/MacroCtrl.md)                                                                                | [finv](../Macros/MacroFinv.md)                                       |
| [elem](../Macros/MacroElem.md)                                                                                | [flo2](../Macros/MacroFlo2.md)                                       |
| [init](../Macros/MacroInit.md) or [pres](../Macros/MacroPres.md)                                             | [flxo](../Macros/MacroFlxo.md) or [flxz](../Macros/MacroFlxz.md)    |
| [perm](../Macros/MacroPerm.md)                                                                                | [fper](../Macros/MacroFper.md)                                       |
| [rlp](../Macros/MacroRlp.md)                                                                                  | [gdpm](../Macros/MacroGdpm.md)                                       |
| [rock](../Macros/MacroRock.md)                                                                                | [hflx](../Macros/MacroHflx.md)                                       |
| [sol](../Macros/MacroSol.md)                                                                                  | [iter](../Macros/MacroIter.md)                                       |
| [time](../Macros/MacroTime.md)                                                                                | [node](../Macros/MacroNode.md) or [nod2](../Macros/MacroNod2.md)    |
| [stop](../Macros/MacroStop.md)                                                                                | [ppor](../Macros/MacroPpor.md)                                       |
|                                                                                                                | [renu](../Macros/MacroRenu.md)                                       |
| dual (* only)                                                                                                  | [rflx](../Macros/MacroRflx.md)                                       |
| dpdp (** only)                                                                                                 | [rxn](../Macros/MacroRxn.md)                                         |
|                                                                                                                | text or comments (#)                                                  |
|                                                                                                                | [trac](../Macros/MacroTrac.md)                                       |
|                                                                                                                | [user](../Macros/MacroUser.md) or userc  |
|                                                                                                                | [vcon](../Macros/MacroVcon.md)                                       |
|                                                                                                                | velo                                       |
|                                                                                                                | [zone](../Macros/MacroZone.md) or [zonn](../Macros/MacroZonn.md)    |


## Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability 


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| title                                                                    | [adif](../Macros/MacroAdif.md)                                       |
| [boun](../Macros/MacroBoun.md) or [flow](../Macros/MacroFlow.md) or [hflx](../Macros/MacroHflx.md) | [cden](../Macros/MacroCden.md)                                       |
| [cond](../Macros/MacroCond.md)                                                                       | [cont](../Macros/MacroCont.md)                                       |
| [coor](../Macros/MacroCoor.md)                                                                       | [eos](../Macros/MacroEos.md)                                         |
| [ctrl](../Macros/MacroCtrl.md)                                                                       | [finv](../Macros/MacroFinv.md)                                       |
| [elem](../Macros/MacroElem.md)                                                                       | [flo2](../Macros/MacroFlo2.md)                                       |
| [init](../Macros/MacroInit.md) or [pres](../Macros/MacroPres.md)                                    | [flxo](../Macros/MacroFlxo.md)                                       |
| [ngas](../Macros/MacroNgas.md)                                                                       | [fper](../Macros/MacroFper.md)                                       |
| [perm](../Macros/MacroPerm.md)                                                                       | [gdpm](../Macros/MacroGdpm.md)                                       |
| [rlp](../Macros/MacroRlp.md)                                                                         | [iter](../Macros/MacroIter.md)                                       |
| [rock](../Macros/MacroRock.md)                                                                       | [node](../Macros/MacroNode.md) or [nod2](../Macros/MacroNod2.md)    |
| [sol](../Macros/MacroSol.md)                                                                         | [ppor](../Macros/MacroPpor.md)                                       |
| [time](../Macros/MacroTime.md)                                                                       | [renu](../Macros/MacroRenu.md)                                       |
| [stop](../Macros/MacroStop.md)                                                                       | [rflx](../Macros/MacroRflx.md)                                       |
|                                                                                                       | [rxn](../Macros/MacroRxn.md)                                         |
| dual (*only)                                                                                          | [szna](../Macros/MacroSzna.md)                                       |
| dpdp (**only)                                                                                         | text or comments (#)                                                  |
|                                                                                                       | [trac](../Macros/MacroTrac.md)                                       |
|                                                                                                       | [user](../Macros/MacroUser.md) or userc  |
|                                                                                                       | [vapl](../Macros/MacroVapl.md)                                       |
|                                                                                                       | [vcon](../Macros/MacroVcon.md)                                       |
|                                                                                                       | velo                                       |
|                                                                                                       | [zone](../Macros/MacroZone.md) or [zonn](../Macros/MacroZone.md)    |



## Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability


| Required Macros |   Optional Macros   | 
|:---------------|:-----------------------|
| title                                    | [bous](../Macros/MacroBous.md)                                       |
| [airwater](../Macros/MacroAirwater.md)                              | [cont](../Macros/MacroCont.md)                                       |
| [boun](../Macros/MacroBoun.md) or [flow](../Macros/MacroFlow.md)   | [eos](../Macros/MacroEos.md)                                         |
| [coor](../Macros/MacroCoor.md)                                      | [exri](../Macros/MacroExri.md)                                       |
| [ctrl](../Macros/MacroCtrl.md)                                      | [finv](../Macros/MacroFinv.md)                                       |
| [elem](../Macros/MacroElem.md)                                      | [flo2](../Macros/MacroFlo2.md)                                       |
| [init](../Macros/MacroInit.md) or [pres](../Macros/MacroPres.md)   | [flxo](../Macros/MacroFlxo.md)                                       |
| [node](../Macros/MacroNode.md) or [nod2](../Macros/MacroNod2.md)   | [fper](../Macros/MacroFper.md)                                       |
| [perm](../Macros/MacroPerm.md)                                      | [gdpm](../Macros/MacroGdpm.md)                                       |
| [rock](../Macros/MacroRock.md)                                      | [head](../Macros/MacroHead.md)                                       |
| [sol](../Macros/MacroSol.md)                                        | [iter](../Macros/MacroIter.md)                                       |
| [time](../Macros/MacroTime.md)                                      | [ppor](../Macros/MacroPpor.md)                                       |
| [stop](../Macros/MacroStop.md)                                      | [pres](../Macros/MacroPres.md)                                       |
|                                                                      | [renu](../Macros/MacroRenu.md)                                       |
| dual (*only)                                                         | [rlp](../Macros/MacroRlp.md)                                         |
| dpdp (only)                                                          | [rxn](../Macros/MacroRxn.md)                                         |
|                                                                      | text or comments                                                      |
|                                                                      | [trac](../Macros/MacroTrac.md)                                       |
|                                                                      | [user](../Macros/MacroUser.md) or userc  |
|                                                                      | [vapl](../Macros/MacroVapl.md)                                       |
|                                                                      | velo                                       |
|                                                                      | [zone](../Macros/MacroZone.md) or [zonn](../Macros/MacroZonn.md)    |



### Required Macros

   * [anpe](../Macros/MacroAnpe.md)
   * [boun](../Macros/MacroBoun.md)
   * [cond](../Macros/MacroCond.md)
   * [coor](../Macros/MacroCoor.md)
   * [elem](../Macros/MacroElem.md)
   * [fdm](../Macros/MacroFdm.md)
   * [flow](../Macros/MacroFlow.md)
   * [ftsc](../Macros/MacroFtsc.md)
   * [hcon](../Macros/MacroHcon.md)
   * [hyco](../Macros/MacroHyco.md)
   * [init](../Macros/MacroInit.md)
   * [ivfc](../Macros/MacroIvfc.md)
   * [perm](../Macros/MacroPerm.md)
   * [pres](../Macros/MacroPres.md)
   * [rock](../Macros/MacroRock.md)
   * [stop](../Macros/MacroStop.md)
   * [time](../Macros/MacroTime.md)



### Optional Macros


   * [adif](../Macros/MacroAdif.md)
   * [airwater or air](../Macros/MacroAirwater.md)
   * [alta](../Macros/MacroAlti.md) (Deprecated)
   * [bous](../Macros/MacroBous.md)
   * [carb](../Macros/MacroCarb.md)
   * [cden](../Macros/MacroCden.md)
   * [cflx](../Macros/MacroCflx.md)
   * [cgdp](../Macros/MacroCgdp.md)
   * [chea](../Macros/MacroChea.md)
   * [conn](../Macros/MacroConn.md)
   * [cont](../Macros/MacroCont.md)
   * [conv](../Macros/MacroConv.md)
   * [ctrl](../Macros/MacroCtrl.md)
   * [dpdp](../Macros/MacroDpdp.md)
   * [dual](../Macros/MacroDual.md)
   * [dvel](../Macros/MacroDvel.md)
   * [eos](../Macros/MacroEos.md)
   * [evap](../Macros/MacroEvap.md)
   * [exrl](../Macros/MacroExrl.md)
   * [finv](../Macros/MacroFinv.md)
   * [flgh](../Macros/MacroFlgh.md)
   * [flo2](../Macros/MacroFlo2.md)
   * [flo3](../Macros/MacroFlo3.md)
   * [floa](../Macros/MacroFloa.md)
   * [flwt](../Macros/MacroFlwt.md)
   * [flxn](../Macros/MacroFlxn.md)
   * [flxo](../Macros/MacroFlxo.md)
   * [flxz](../Macros/MacroFlxz.md)
   * [fper](../Macros/MacroFper.md)
   * [frlp](../Macros/MacroFrlp.md)
   * [gdkm](../Macros/MacroGdkm.md)
   * [gdpm](../Macros/MacroGdpm.md)
   * [grad](../Macros/MacroGrad.md)
   * [head](../Macros/MacroHead.md)
   * [hflx](../Macros/MacroHflx.md)
   * [hist](../Macros/MacroHist.md)
   * [ice](../Macros/MacroIce.md) or [meth](../Macros/MacroIce.md)
   * [intg](../Macros/MacroIntg.md)
   * [imex](../Macros/MacroImex.md)
   * [impf](../Macros/MacroImpf.md)
   * [isot](../Macros/MacroIsot.md)
   * [iter](../Macros/MacroIter.md)
   * [itfc](../Macros/MacroItfc.md)
   * [ittm](../Macros/MacroIttm.md)
   * [itup](../Macros/MacroItup.md)
   * [iupk](../Macros/MacroIupk.md)
   * [mdnode](../Macros/MacroMdnode.md)
   * [mptr](../Macros/MacroMptr.md)
   * [nfinv](../Macros/MacroNfinv.md)
   * [ngas](../Macros/MacroNgas.md)
   * [nobr](../Macros/MacroNobr.md)
   * [nod2](../Macros/MacroNod2.md)
   * [nod3](../Macros/MacroNod3.md)
   * [node](../Macros/MacroNode.md)
   * [nrst](../Macros/MacroNrst.md)
   * [para](../Macros/MacroPara.md)
   * [pest](../Macros/MacroPest.md)
   * [phys](../Macros/MacroPhys.md)
   * [ppor](../Macros/MacroPpor.md)
   * [ptrk](../Macros/MacroPtrk.md)
   * [renu](../Macros/MacroRenu.md)
   * [rest](../Macros/MacroRest.md)
   * [rflo](../Macros/MacroRflo.md)
   * [rflx](../Macros/MacroRflx.md) (Deprecated)
   * [rich](../Macros/MacroRich.md)
   * [rive](../Macros/MacroRive.md) or [well](../Macros/MacroRive.md)
   * [rlp](../Macros/MacroRlp.md)
   * [rlpm](../Macros/MacroRlpm.md)
   * [rxn](../Macros/MacroRxn.md)
   * [sol](../Macros/MacroSol.md)
   * [sptr](../Macros/MacroSptr.md)
   * [stea](../Macros/MacroStea.md)
   * [strs](../Macros/MacroStrs.md)
   * [subm](../Macros/MacroSubm.md)
   * [svar](../Macros/MacroSvar.md)
   * [szna](../Macros/MacroSzna.md)
   * [text](../Macros/MacroText.md)
   * [thic](../Macros/MacroThic.md)
   * [trac](../Macros/MacroTrac.md)
   * [trxn](../Macros/MacroTrxn.md)
   * [user](../Macros/MacroUser.md)
   * [vapl](../Macros/MacroVapl.md)
   * [vbou](../Macros/MacroVbou.md)
   * [vcon](../Macros/MacroVcon.md)
   * [weli](../Macros/MacroWeli.md)
   * [wflo](../Macros/MacroWflo.md)
   * [wgtu](../Macros/MacroWgtu.md)
   * [wtsi](../Macros/MacroWtsi.md)
   * [zeol](../Macros/MacroZeol.md)
   * [zneg](../Macros/MacroZneg.md)
   * [zone](../Macros/MacroZone.md)
   * [zonn](../Macros/MacroZonn.md)
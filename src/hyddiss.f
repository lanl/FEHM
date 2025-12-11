      subroutine hyddiss(ndummy,iz)
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To compute the dissociation models for hydrate.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation:  8-05-2005,  Programmer: George Zyvoloski
CD2
CD2 $Log:   /pvcs.config/fehm90/src/hyddiss.f_a  $
CD2                                     
C**********************************************************************
c
c calculates dissociation rates and derivatives
c
      use comai
      use combi, only : sx1
      use comdi, only : ps
      use commeth
      use comdti
      use comki
      implicit none

      integer iz,ndummy,i,idissd,mi,ieosd,it,ir,j,num_models,ireg

      logical null1,ex

      real*8 var1,var2,var3,var4
      real*8 prop,der1,der2,der3,der4
      real*8 tp,pt,dtp,dpt2,gvar,dgterm3,dgterm4,actterm,dactterm2
      real*8 aterm,hterm,dhterm4,pdif,pterm,dpterm1,dpterm2
      real*8 porterm,mol_hyd,wterm,dwterm3,gvar1,gterm,u_tol
      real*8 gas_const,temp_conv,v_tol,vard,eact2,den_hyd,den_hyd0
      real*8 dis_tol,prop_tol
      real*8 c_eq_h_lw,c_eq_lw_v,flux_kiyono_1,flux_kiyono
      parameter (gas_const = 8.314d0, temp_conv = 273.15d0)
      parameter (v_tol = 1.d-6,u_tol = 1.d0)
      parameter (mol_hyd = 0.1195d0, den_hyd0 = 7866.d0)
      parameter (dis_tol=0.02,  prop_tol =1.0d-40)
      
c     da - grain area term
c     oe - exponent of grain area term
c     pe - exponent of fugacity difference
c     le - exponent of hydrate saturation term
c     me - exponent of water saturation term
c     ne - exponent of gas saturation term
c     e_act - activation energy  (two terms)
c     afhydh - activity fequency term 

c     dar - grain area term (reformation)
c     oer - exponent of grain area term (reformation)
c     per - exponent of fugacity difference (reformation)
c     ler - exponent of hydrate saturation term (reformation)
c     mer - exponent of water saturation term (reformation)
c     ner - exponent of gas saturation term (reformation)
c     e_actr - activation energy (reformation or growth)
c     afhydr - activity fequency term (reformation or growth)

      if(iz.eq.0) then
c     first allocate memory for models      
         allocate (idissp(n0))
         allocate (idisst(max_hyd_model))

         allocate (diss1f(max_hyd_model))
         allocate (diss2f(max_hyd_model))
         allocate (diss3f(max_hyd_model))
         allocate (diss4f(max_hyd_model))
         allocate (diss5f(max_hyd_model))
         allocate (diss6f(max_hyd_model))
         allocate (diss7f(max_hyd_model))
         allocate (diss8f(max_hyd_model))
         allocate (diss9f(max_hyd_model))
         allocate (form1f(max_hyd_model))
         allocate (form2f(max_hyd_model))
         allocate (form3f(max_hyd_model))
         allocate (form4f(max_hyd_model))
         allocate (form5f(max_hyd_model))
         allocate (form6f(max_hyd_model))
         allocate (form7f(max_hyd_model))
         allocate (form8f(max_hyd_model))
         allocate (form9f(max_hyd_model))
         
c     read in data

c     check for read from other file

         i = 0
         j = 0
         ex = .false.
 10      continue
         read(inpt,'(a80)') wdd1
         if(.not. null1(wdd1)) then
            backspace inpt
            read(inpt,*) idissd
            backspace inpt
            i = i+1
            if (i.gt.max_hyd_model) then 
c     error condition  to many models
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,*) 'error ',
     &              'in hydiss: too many models (line 111 hydiss.f)'
               if (iatty.gt.0) write(iatty,*) 'error in hydiss: ',
     &              'too many models (line 111 hydiss.f)'
            endif
            if (idissd .eq. 1) then
c     this is previous model4, the simplest rate equation  (temperature based) 
c     diss1f(i) = da - grain area term
               read(inpt,*) idisst(i),diss1f(i)
               
            else  if (idissd .eq. 2) then
c     this is the kim-bishnoi model and is used for disociation 
c     there is no forming  
c     diss1f(i) da - grain area term
c     diss7f(i) e_act - activation energy 
c     diss8f(i) afhydh - activity fequency term
               read(inpt,*) idisst(i),diss1f(i),diss7f(i),diss8f(i)
               
            else  if (idissd .eq. 3) then
c     this is the kim-bishnoi model (dissociation)
c     and general formation model 
c     dissociation:
c     diss1f(i) da - grain area term
c     diss7f(i) e_act - activation energy 
c     diss8f(i) afhydh - activity fequency term
c     forming:
c     form1f(i) dar - grain area term
c     form2f(i) oer - exponent of grain area term
c     form3f(i) per - exponent of fugacity difference
c     form4f(i) ler - exponent of hydrate saturation term
c     form5f(i) mer - exponent of water saturation term
c     form6f(i) ner - exponent of gas saturation term
c     form7f(i) e_actr - activation energy(1)
c     form8f(i) e_actr - activation energy(2)  
c     form9f(i) afhydhr - activity fequency term
               read(inpt,*) idisst(i),diss1f(i),diss7f(i),diss8f(i),
     &              form1f(i),form2f(i),form3f(i),form4f(i),
     &              form5f(i), form6f(i), form7f(i), form8f(i),
     &              form9f(i)
               if(form2f(i).lt.0.0) then
                  form2f(i)=abs(form2f(i))
                  form1f(i) = 1.0/form1f(i)
               endif  
            else  if (idissd .eq. 4) then
c     c Kiyono type hydrate dissociation and 
c     general model (sakamoto) for hydrate forming 
c     same input as model 3 
               read(inpt,*) idisst(i),diss1f(i),diss7f(i),diss8f(i),
     &              form1f(i),form2f(i),form3f(i),form4f(i),
     &              form5f(i), form6f(i), form7f(i), form8f(i),
     &              form9f(i)
               if(form2f(i).lt.0.0) then
                  form2f(i)=abs(form2f(i))
                  form1f(i) = 1.0/form1f(i)
               endif  
            else  if (idissd .eq. 5) then
c     this is the kim-bishnoi model (dissociation)
c     and general formation model 
c     dissociation:
c     diss1f(i) da - grain area term
c     diss7f(i) e_act - activation energy 
c     diss8f(i) afhydh - activity fequency term
c     forming:
c     form1f(i) dar - grain area term
c     form2f(i) oer - exponent of grain area term
c     form3f(i) per - exponent of fugacity difference
c     form4f(i) ler - exponent of hydrate saturation term
c     form5f(i) mer - exponent of water saturation term
c     form6f(i) ner - exponent of gas saturation term
c     form7f(i) e_actr - activation energy(1)
c     form8f(i) e_actr - activation energy(2)  
c     form9f(i) afhydhr - activity fequency term
               read(inpt,*) idisst(i),diss1f(i),diss7f(i),diss8f(i),
     &              form1f(i),form2f(i),form3f(i),form4f(i),
     &              form5f(i), form6f(i), form7f(i), form8f(i),
     &              form9f(i)
               if(form2f(i).lt.0.0) then
                  form2f(i)=abs(form2f(i))
                  form1f(i) = 1.0/form1f(i)
               endif
            else  if (idissd .eq. 6) then
c     this is the general formation model  for dissociation and forming 
c     dissociation:
c     diss1f(i) da - grain area term
c     diss2f(i) oe - exponent of grain area term
c     diss3f(i) pe - exponent of fugacity difference
c     diss4f(i) le - exponent of hydrate saturation term
c     diss5f(i) me - exponent of water saturation term
c     diss6f(i) ne - exponent of gas saturation term
c     diss7f(i) e_act - activation energy (1)
c     diss8f(i) e_act - activation energy (2)  
c     diss9f(i) afhydh - activity fequency term
c     forming:
c     form1f(i) dar - grain area term
c     form2f(i) oer - exponent of grain area term
c     form3f(i) per - exponent of fugacity difference
c     form4f(i) ler - exponent of hydrate saturation term
c     form5f(i) mer - exponent of water saturation term
c     form6f(i) ner - exponent of gas saturation term
c     form7f(i) e_actr - activation energy (1) 
c     form8f(i) e_actr - activation energy (2) 
c     form9f(i) afhydhr - activity fequency term
               read(inpt,*) idisst(i),diss1f(i),diss2f(i),diss3f(i),
     &              diss4f(i), diss5f(i), diss6f(i), 
     &              diss7f(i), diss8f(i),diss9f(i),
     &              form1f(i),form2f(i),form3f(i),
     &              form4f(i), form5f(i), form6f(i), 
     &              form7f(i), form8f(i),form9f(i)
               
            endif 
            go to 10
         endif                
c     
c     return if read data from a file
c     
         if(ex) return
         
         
c     read in nodal capillary type
c gaz 110518   macroread(7) changed macroread(25)      
         narrays = 1
         itype(1) = 4
         default(1) = 1
         macro = "methdisv "
         igroup = 2
         call initdata2( inpt, ischk, n0, narrays,
     2      itype, default, macroread(25), macro(1:4), igroup, ireturn,
     3      i4_1=idissp(1:n0) )

         macroread(25) = .TRUE.
         
         do i=1,n0
            if(idisst(idissp(i)).lt.1) then
c     error condition  to many models
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,*) 'error i',
     &              'n hydiss: invalid model number (line 202 hydiss.f)'
               if (iatty.gt.0) write(iatty,*) 'error in hydiss: ',
     &              'invalid model number (line 202 hydiss.f)'       
            else if(idisst(idissp(i)).gt.5) then
c     error condition  
            endif
         enddo
         
         
      else
c     
c     calculate  hydrate rates and derivatives
c     
         do i = 1,neq
            mi = i+ndummy
c     it =  model number
            it = idissp(mi)
c     idissd =  model type            
            if(it.eq.0) then
               idissd = 0
            else
               idissd = idisst(it)
            endif
            if(idissd.eq.1) then
c     
c     simple temperature based rate equation (dissociation and forming)
c     diss1f(it) da - grain area term
c     
               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               da = diss1f(it)
               volh = sx1(mi)
               porh = ps(mi) 
               tp = (log(var1)-ptc2)/ptc1 - temp_conv
               dtp = 1.0/(ptc1*var1)
               me =(var2-tp)
               if(me.lt.0) then
c     growing (need gas and water)
c     gas fraction + water term
                  gvar = (1.0-var3-var4)*var3
                  dgterm3 =  -1.0*var3 + (1.0-var3-var4)
                  dgterm4 =  -1.0*var3
                  oe = volh*porh*da
                  prop = oe*gvar*me
                  der1 = -oe*gvar*dtp
                  der2 = oe*gvar
                  der3 = oe*me*dgterm3
                  der4 = oe*me*dgterm4
               else
c     dissociation (need hydrate)
                  gvar = var4
c     
                  if (gvar .ge. 1.0) then
                     gvar = 1.0
	          end if
                  if (gvar .le. 0.0) then
                     gvar = 0.0
	          end if
c     
                  dgterm3 =  0.0
                  dgterm4 =  1.0
                  oe = volh*porh*da
                  prop = oe*gvar*me
                  der1 = -oe*gvar*dtp
                  der2 = oe*gvar
                  der3 = oe*me*dgterm3
                  der4 = oe*me*dgterm4              
               end if
               
               skhyd_temp(mi) = prop
               dskhydpf(mi) = der1
               dskhydtf(mi) = der2
               dskhydwf(mi) = der3
               dskhydmf(mi) = der4  
               
            else if(idissd.eq.2) then
c     this is the kim-bishnoi model for hydrate disociation and no hydrate forming 
c     diss1f(it) da - grain area term
c     diss7f(it) e_act - activation energy 
c     diss8f(it) afhydh - activity fequency term
               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               da     = diss1f(it)
               e_act1 = diss7f(it)
               afhyd  = diss8f(it)
               pt = exp((ptc1*(var2+temp_conv))+ptc2)
               if (var1.le.pt) then
                  dpt2 = pt*ptc1
c     hydrate can dissociate 
                  eact2 = e_act1/gas_const
                  den_hyd = 1.0
c     activity term
                  actterm   =  exp(-eact2/(var2+temp_conv))
                  dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
c     area term         
                  aterm  =  da
c     hydrate fraction term          
                  hterm  =  var4
                  dhterm4 = 1.0
c     pressure difference term         
                  pdif   =  pt-var1
                  pterm   =  pdif
                  dpterm1 =  -1.0
                  dpterm2 =  dpterm1*dpt2
c     
c     porosity, perm, density term
                  volh = sx1(mi)
                  porh = ps(mi)
                  porterm = mol_hyd*porh*den_hyd*afhyd*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = 0.0d0
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*pterm*dhterm4
                  
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4
                  
               else
c     no hydrate forming allowed in this model  
                  skhyd_temp(mi) = 0.0d0
                  dskhydpf(mi) = 0.0d0
                  dskhydtf(mi) = 0.0d0
                  dskhydwf(mi) = 0.0d0
                  dskhydmf(mi) = 0.0d0  
               endif
            else if(idissd.eq.3) then
c     this is the kim-bishnoi model for hydrate disociation 
c     general model (sakamoto) for hydrate forming 

               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               porh = ps(mi)
               volh = sx1(mi)
               pt = exp((ptc1*(var2+temp_conv))+ptc2)
               if (var1.le.pt) then
c     diss1f(it) da - grain area term
c     diss7f(it) e_act - activation energy 
c     diss8f(it) afhydh - activity fequency term
                  da     = diss1f(it)
                  e_act1 = diss7f(it)
                  afhyd  = diss8f(it)
                  dpt2 = pt*ptc1
c     hydrate can dissociate 
                  eact2 = e_act1/gas_const
                  den_hyd = 1.0
c     activity term
                  actterm   =  exp(-eact2/(var2+temp_conv))
                  dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
c     area term         
                  aterm  =  da
c     hydrate fraction term          
                  hterm  =  var4
                  dhterm4 = 1.0
c     pressure difference term         
                  pdif   =  pt-var1
                  pterm   =  pdif
                  dpterm1 =  -1.0
                  dpterm2 =  dpterm1*dpt2
c     
c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhyd*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = 0.0d0
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*pterm*dhterm4
                  
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4
                  
               else
c     general model for hydrate forming
c     Sakamoto type 
                  dar 	    = form1f(it)
                  oer 	    = form2f(it)
                  per 	    = form3f(it)
                  ler 	    = form4f(it)
                  mer 	    = form5f(it)
                  ner 	    = form6f(it)
                  e_act1r     = form7f(it)
                  e_act2r     = form8f(it)
                  afhydr      = form9f(it)
                  
                  den_hyd = den_hyd0 

c     activity term
                  actterm   =  exp(e_act1r/(var2+temp_conv)+e_act2r)
                  dactterm2 =  -actterm*(e_act1r/(var2+temp_conv)**2)
                  aterm  =  dar**oer
c     hydrate fraction term (included hydrate frac from start of eq.)
                  if(ler+1.0.gt.0.0) then
                     if(var4.lt.0) then
                        vard = 0.0
                     else
                        vard = var4
                     endif
                     hterm   =  vard**(ler+1.0)
                     if(ler+1.0.eq.1) then
                        dhterm4 = 1.d0
                     else
                        dhterm4 =  (ler+1.0)*vard**(ler)
                     endif
                  else
                     hterm = 1.0
                     dhterm4 = 0.0
                  endif
                  if(mer.gt.0.0) then
                     if(var3.lt.0) then
                        vard = 0.0
                     else
                        vard = var3
                     endif
                     wterm   =  vard**mer
                     if(mer.ne.1) then
                        dwterm3 =  mer*vard**(mer-1.0)
                     else
                        dwterm3 = 1.0d0
                     endif
                  else
                     wterm = 1.0
                     dwterm3 = 0.0
                  endif
c     gas fraction term
                  gvar = 1.0-var3-var4
                  if(ner.gt.0.0) then
                     if(gvar.lt.0) then
                        gvar1 = 0.0
                     else
                        gvar1 = gvar
                     endif
                     gterm = gvar1**ner			
                     if(ner.ne.1) then
                        dgterm3 =  -ner*gvar1**(ner-1.0)
                        dgterm4 =  dgterm3
                     else
                        dgterm3 = 1.0d0
                        dgterm4 = 1.0d0
                     endif
                  else
                     gterm = 1.0
                     dgterm3 = 0.0
                     dgterm4 = 0.0
                  endif
c     pressure (fugacity) term
c     pdif need to be able to change sign
c     might restrict pe to greater than or equal to 1
                  pdif   =  pt-var1
                  pterm   =  pdif**per
                  dpterm1 =  -per*pdif**(per-1.0)
                  dpterm2 =  -dpterm1*dpt2

c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhydr*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*wterm*gterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &                 + wterm*dgterm3)
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &                 + hterm*dgterm4)   
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4  
               endif

            else if(idissd.eq.4) then
c     Kiyono type hydrate dissociation and 
c     general model (sakamoto) for hydrate forming 
               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               porh = ps(mi)
               volh = sx1(mi)
               pt = exp((ptc1*(var2+temp_conv))+ptc2)
               if (var1.le.pt) then
c     Kiyono type hydrate dissociation
c     c  tenma 2005/Jan/18
c     Kiyono type  ( diss_type is not 1.0 )

                  den_hyd = den_hyd0

c     porosity, perm, density term
                  porterm = mol_hyd*porh*volh*u_tol
c     Flux calculation  < tenma 2005/Jan/18 >
                  c_eq_h_lw=exp(-152.777+7478.8/(var2+temp_conv)+
     &                 20.6794*log(var2+temp_conv)+0.75316*log(10*pt))
c     
                  c_eq_lw_v=exp(-152.777+7478.8/(var2+temp_conv)+
     &                 20.6794*log(var2+temp_conv)+0.75316*log(10*var1))
c     
                  flux_kiyono_1=afhyd*gas_const*log(c_eq_h_lw/c_eq_lw_v)
                  flux_kiyono=flux_kiyono_1*
     &                 exp(e_act1/(var2+temp_conv))*(var2+temp_conv)
c     
c     rate term is made up of individual terms
                  prop = da*porterm*var4*flux_kiyono
c     pressure derivative  (der1)
                  der1 = 0.0
c     temperature derivative  (der2)
                  der2 = da*porterm*var4*flux_kiyono_1*
     &                 (1-(e_act1/(var2+temp_conv)))*
     &                 exp(e_act1/(var2+temp_conv))
c     water fraction derivative  (der3)
                  der3 = 0.0
c     hydrate fraction derivative  (der4)   
                  der4 = da*porterm*flux_kiyono 
               else
c     general model for hydrate forming
c     Sakamoto type 
                  dar 	    = form1f(it)
                  oer 	    = form2f(it)
                  per 	    = form3f(it)
                  ler 	    = form4f(it)
                  mer 	    = form5f(it)
                  ner 	    = form6f(it)
                  e_act1r     = form7f(it)
                  e_act2r     = form8f(it)
                  afhydr      = form9f(it)
                  
                  den_hyd = den_hyd0 

c     activity term
                  actterm   =  exp(e_act1r/(var2+temp_conv)+e_act2r)
                  dactterm2 =  -actterm*(e_act1r/(var2+temp_conv)**2)
                  aterm  =  dar**oer
c     hydrate fraction term (included hydrate frac from start of eq.)
                  if(ler+1.0.gt.0.0) then
                     if(var4.lt.0) then
                        vard = 0.0
                     else
                        vard = var4
                     endif
                     hterm   =  vard**(ler+1.0)
                     if(ler+1.0.eq.1) then
                        dhterm4 = 1.d0
                     else
                        dhterm4 =  (ler+1.0)*vard**(ler)
                     endif
                  else
                     hterm = 1.0
                     dhterm4 = 0.0
                  endif
                  if(mer.gt.0.0) then
                     if(var3.lt.0) then
                        vard = 0.0
                     else
                        vard = var3
                     endif
                     wterm   =  vard**mer
                     if(mer.ne.1) then
                        dwterm3 =  mer*vard**(mer-1.0)
                     else
                        dwterm3 = 1.0d0
                     endif
                  else
                     wterm = 1.0
                     dwterm3 = 0.0
                  endif
c     gas fraction term
                  gvar = 1.0-var3-var4
                  if(ner.gt.0.0) then
                     if(gvar.lt.0) then
                        gvar1 = 0.0
                     else
                        gvar1 = gvar
                     endif
                     gterm = gvar1**ner			
                     if(ner.ne.1) then
                        dgterm3 =  -ner*gvar1**(ner-1.0)
                        dgterm4 =  dgterm3
                     else
                        dgterm3 = 1.0d0
                        dgterm4 = 1.0d0
                     endif
                  else
                     gterm = 1.0
                     dgterm3 = 0.0
                     dgterm4 = 0.0
                  endif
c     pressure (fugacity) term
c     pdif need to be able to change sign
c     might restrict pe to greater than or equal to 1
                  pdif   =  pt-var1
                  pterm   =  pdif**per
                  dpterm1 =  -per*pdif**(per-1.0)
                  dpterm2 =  -dpterm1*dpt2

c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhydr*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*wterm*gterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &                 + wterm*dgterm3)
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &                 + hterm*dgterm4)   
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4  
               endif

            else if(idissd.eq.5) then
c     this is the kim-bishnoi model for hydrate disociation 
c     general model (sakamoto) for hydrate forming 

               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               porh = ps(mi)
               volh = sx1(mi)
               pt = exp((ptc1*(var2+temp_conv))+ptc2)
               if (var1.le.pt) then
c     diss1f(it) da - grain area term
c     diss7f(it) e_act - activation energy 
c     diss8f(it) afhydh - activity fequency term
                  da     = diss1f(it)
                  e_act1 = diss7f(it)
                  afhyd  = diss8f(it)
                  dpt2 = pt*ptc1
c     hydrate can dissociate 
                  eact2 = e_act1/gas_const
                  den_hyd = 1.0
c     activity term
                  actterm   =  exp(-eact2/(var2+temp_conv))
                  dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
c     area term         
                  aterm  =  da
c     hydrate fraction term => tenma modify 2005/11/07          
                  vard  =  var4
                  if (vard.gt.frachyd_1(i)+dis_tol) then
                     hterm = vard - frachyd_1(i)
                     dhterm4 = 1.0
                  else
c     gaz 7-17-2006 (note that this produces a discontinuous function)
                     hterm = vard
                     dhterm4 = 1.0
c     hterm = 0.0
c     dhterm4 = 0.0
                  end if
c     
c     pressure difference term         
                  pdif   =  pt-var1
                  pterm   =  pdif
                  dpterm1 =  -1.0
                  dpterm2 =  dpterm1*dpt2
c     
c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhyd*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = 0.0d0
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*pterm*dhterm4
                  if(abs(prop).ge.prop_tol) then
                     skhyd_temp(mi) = prop
                     dskhydpf(mi) = der1
                     dskhydtf(mi) = der2
                     dskhydwf(mi) = der3
                     dskhydmf(mi) = der4
                  else
                     skhyd_temp(mi) = 0.0
                     dskhydpf(mi) = 0.0
                     dskhydtf(mi) = 0.0
                     dskhydwf(mi) = 0.0
                     dskhydmf(mi) = 0.0
                  endif
                  
               else
c     general model for hydrate forming
c     Sakamoto type 
                  dar 	    = form1f(it)
                  oer 	    = form2f(it)
                  per 	    = form3f(it)
                  ler 	    = form4f(it)
                  mer 	    = form5f(it)
                  ner 	    = form6f(it)
                  e_act1r     = form7f(it)
                  e_act2r     = form8f(it)
                  afhydr      = form9f(it)
                  
                  den_hyd = den_hyd0 

c     activity term
                  actterm   =  exp(e_act1r/(var2+temp_conv)+e_act2r)
                  dactterm2 =  -actterm*(e_act1r/(var2+temp_conv)**2)
                  aterm  =  dar**oer
c     hydrate fraction term (included hydrate frac from start of eq.)
                  if(ler+1.0.gt.0.0) then
                     if(var4.lt.0) then
                        vard = 0.0
                     else
                        vard = var4
                     endif
                     hterm   =  vard**(ler+1.0)
                     if(ler+1.0.eq.1) then
                        dhterm4 = 1.d0
                     else
                        dhterm4 =  (ler+1.0)*vard**(ler)
                     endif
                  else
                     hterm = 1.0
                     dhterm4 = 0.0
                  endif
                  if(mer.gt.0.0) then
                     if(var3.lt.0) then
                        vard = 0.0
                     else
                        vard = var3
                     endif
                     wterm   =  vard**mer
                     if(mer.ne.1) then
                        dwterm3 =  mer*vard**(mer-1.0)
                     else
                        dwterm3 = 1.0d0
                     endif
                  else
                     wterm = 1.0
                     dwterm3 = 0.0
                  endif
c     gas fraction term
                  gvar = 1.0-var3-var4
                  if(ner.gt.0.0) then
                     if(gvar.lt.0) then
                        gvar1 = 0.0
                     else
                        gvar1 = gvar
                     endif
                     gterm = gvar1**ner			
                     if(ner.ne.1) then
                        dgterm3 =  -ner*gvar1**(ner-1.0)
                        dgterm4 =  dgterm3
                     else
                        dgterm3 = 1.0d0
                        dgterm4 = 1.0d0
                     endif
                  else
                     gterm = 1.0
                     dgterm3 = 0.0
                     dgterm4 = 0.0
                  endif
c     pressure (fugacity) term
c     pdif need to be able to change sign
c     might restrict pe to greater than or equal to 1
                  pdif   =  pt-var1
                  pterm   =  pdif**per
                  dpterm1 =  -per*pdif**(per-1.0)
                  dpterm2 =  -dpterm1*dpt2

c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhydr*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*wterm*gterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &                 + wterm*dgterm3)
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &                 + hterm*dgterm4)
                  if(abs(prop).ge.prop_tol) then   
                     skhyd_temp(mi) = prop
                     dskhydpf(mi) = der1
                     dskhydtf(mi) = der2
                     dskhydwf(mi) = der3
                     dskhydmf(mi) = der4  
                  else
                     skhyd_temp(mi) = 0.0
                     dskhydpf(mi) = 0.0
                     dskhydtf(mi) = 0.0
                     dskhydwf(mi) = 0.0
                     dskhydmf(mi) = 0.0 
                  endif
                  
               endif

            else if(idissd.eq.6) then
c     this is the kim-bishnoi model for hydrate disociation 
c     general model (sakamoto) for hydrate forming 

               var1 = phimeth(mi)
               var2 = tmeth(mi)
               var3 = fracw(mi)
               var4 = frachyd(mi)
               porh = ps(mi)
               volh = sx1(mi)
               pt = exp((ptc1*(var2+temp_conv))+ptc2)
               if (var1.le.pt) then
c     diss1f(it) da - grain area term
c     diss7f(it) e_act - activation energy 
c     diss8f(it) afhydh - activity fequency term
                  da     = diss1f(it)
                  e_act1 = diss7f(it)
                  afhyd  = diss8f(it)
                  dpt2 = pt*ptc1
c     hydrate can dissociate 
                  eact2 = e_act1/gas_const
                  den_hyd = 1.0
c     activity term
                  actterm   =  exp(-eact2/(var2+temp_conv))
                  dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
c     area term         
                  aterm  =  da
c     hydrate fraction term => tenma modify 2005/11/07          
                  vard  =  var4
                  if (vard.gt.frachyd_1(i)+dis_tol) then
                     hterm = vard - frachyd_1(i)
                  end if
                  if (vard.gt.frachyd_1(i) 
     &                 .and. vard.le.frachyd_1(i)+dis_tol) then
                     hterm = dis_tol
                  end if
                  if (vard.le.frachyd_1(i) ) then
                     hterm = vard
                  end if
                  dhterm4 = 1.0
c     pressure difference term         
                  pdif   =  pt-var1
                  pterm   =  pdif
                  dpterm1 =  -1.0
                  dpterm2 =  dpterm1*dpt2
c     
c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhyd*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = 0.0d0
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*pterm*dhterm4
                  
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4
                  
               else
c     general model for hydrate forming
c     Sakamoto type 
                  dar 	    = form1f(it)
                  oer 	    = form2f(it)
                  per 	    = form3f(it)
                  ler 	    = form4f(it)
                  mer 	    = form5f(it)
                  ner 	    = form6f(it)
                  e_act1r     = form7f(it)
                  e_act2r     = form8f(it)
                  afhydr      = form9f(it)
                  
                  den_hyd = den_hyd0 

c     activity term
                  actterm   =  exp(e_act1r/(var2+temp_conv)+e_act2r)
                  dactterm2 =  -actterm*(e_act1r/(var2+temp_conv)**2)
                  aterm  =  dar**oer
c     hydrate fraction term (included hydrate frac from start of eq.)
                  if(ler+1.0.gt.0.0) then
                     if(var4.lt.0) then
                        vard = 0.0
                     else
                        vard = var4
                     endif
                     hterm   =  vard**(ler+1.0)
                     if(ler+1.0.eq.1) then
                        dhterm4 = 1.d0
                     else
                        dhterm4 =  (ler+1.0)*vard**(ler)
                     endif
                  else
                     hterm = 1.0
                     dhterm4 = 0.0
                  endif
                  if(mer.gt.0.0) then
                     if(var3.lt.0) then
                        vard = 0.0
                     else
                        vard = var3
                     endif
                     wterm   =  vard**mer
                     if(mer.ne.1) then
                        dwterm3 =  mer*vard**(mer-1.0)
                     else
                        dwterm3 = 1.0d0
                     endif
                  else
                     wterm = 1.0
                     dwterm3 = 0.0
                  endif
c     gas fraction term
                  gvar = 1.0-var3-var4
                  if(ner.gt.0.0) then
                     if(gvar.lt.0) then
                        gvar1 = 0.0
                     else
                        gvar1 = gvar
                     endif
                     gterm = gvar1**ner			
                     if(ner.ne.1) then
                        dgterm3 =  -ner*gvar1**(ner-1.0)
                        dgterm4 =  dgterm3
                     else
                        dgterm3 = 1.0d0
                        dgterm4 = 1.0d0
                     endif
                  else
                     gterm = 1.0
                     dgterm3 = 0.0
                     dgterm4 = 0.0
                  endif
c     pressure (fugacity) term
c     pdif need to be able to change sign
c     might restrict pe to greater than or equal to 1
                  pdif   =  pt-var1
                  pterm   =  pdif**per
                  dpterm1 =  -per*pdif**(per-1.0)
                  dpterm2 =  -dpterm1*dpt2

c     porosity, perm, density term
                  porterm = mol_hyd*porh*den_hyd*afhydr*volh*aterm

c     rate term is made up of individual terms
                  prop = porterm*actterm*hterm*wterm*gterm*pterm
c     pressure derivative  (der1)
                  der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c     temperature derivative  (der2)
                  der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &                 + dactterm2*pterm)
c     water fraction derivative  (der3)
                  der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &                 + wterm*dgterm3)
c     hydrate fraction derivative  (der4)   
                  der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &                 + hterm*dgterm4)   
                  skhyd_temp(mi) = prop
                  dskhydpf(mi) = der1
                  dskhydtf(mi) = der2
                  dskhydwf(mi) = der3
                  dskhydmf(mi) = der4  

               endif

            else if(idissd.eq.7) then
c     general hydrate  model for hydrate dissociation and 
c     general hydrate forming model
            endif

         enddo

      endif

      
      return
      end

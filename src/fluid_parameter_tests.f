	module comHT_test
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************


	real*8, allocatable :: tlist(:)
      real*8, allocatable :: ntHT2(:)
      real*8, allocatable :: slist(:)
	real*8, allocatable :: plist(:)
      real*8, allocatable :: pclist(:)
      real*8 tol_var
      real*8 denr_simple,cpr_simple,ps_simple,vol_simple 
      real*8 dmef_0,dmpf_0,dmc_0,depf_0,deef_0,dec_0,dcp_0,dce_0
      real*8 dcc_0, dtpa_0, dtpac_0, dtpae_0, a_1,a_2,a_3,dstm_0 
      real*8 dil_0,dilp_0,dile_0,enlf_0,dglp_0,dgle_0,delf_0,delef_0
      real*8 delcf_0,dclf_0,dclef_0,dclcf_0,dilc_0,dglc_0
      real*8 div_0,divp_0,dive_0,envf_0,dgvp_0,dgve_0,devf_0,devef_0
      real*8 devcf_0,dcvf_0,dcvef_0,dcvcf_0,divc_0,dgvc_0
      real*8 rolf_0,rovf_0,cnlf_0,cnvf_0,enva_0
      real*8 denvap_0, denvae_0, denvac_0
      real*8 var_h2o(25),var_gas(25),var_ngas(25),var_awh_param(30)
      parameter (denr_simple = 2500.d0, cpr_simple = 1000.d0)
      parameter (ps_simple = 1.0d0, vol_simple = 1.d0)
      parameter(tol_var = 1.d-20)
      integer nHTlist, npHTlist, ntHTlist, npcHTlist, nsHTlist 
      integer inctrl, intest ,ierrHT
      real*8 pl_wrt, pcl_wrt, sl_wrt, tl_wrt, t_last, s_last
	end module comHT_test  
    
      subroutine fluid_parameter_tests(iflg) 
c gaz 190725
c this subroutine is used to write out fluid properties
c
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comji
      use davidi
      use comflow
      use commass_AWH
      use comHT_test

      implicit none
c
      integer neqp1,i,j,k,iflg,neq_save,i_awh1,i_awh2
      integer open_file
      real*8 sx1d,tmche,tolerance,psatl,pv,dpsatt,dpsats
      real*8 phi_save,pci_save,t_save,s_save
      if(n_awh_test.eq.0) return
      if(iflg.eq.0) then
c read input (phi,s,t,pci)
       read(inpt,*) n_awh_test
       if(.not.allocated(awh_var_test)) then
        allocate(awh_var_test(n_awh_test,4))
        do i = 1, n_awh_test 
         read(inpt,*) (awh_var_test(i,j), j = 1,4)
        enddo
        continue
       else
        do i = 1, n_awh_test 
         read(inpt,*) (awh_var_test(i,j), j = 1,4)
        enddo
        continue
      endif
      elseif(iflg.eq.1) then
       do i = 1, n_awh_test 
        phi(i) =  awh_var_test(i,1)
        s(i) =  awh_var_test(i,2)
        t(i) =  awh_var_test(i,3)
        pci(i) =  awh_var_test(i,4)
        if(s(i).le.0.d0) then
         ieos(i) = 3
        else if(s(i).ge.1.d0) then
         ieos(i) = 1
        else 
         ieos(i) = 2
          pv = phi(i) - pci(i)
          t(i)= psatl(pv,0.0d0,dpcef(i),dpsatt,dpsats,1,an(i))   
        endif                  
       enddo
       continue
      elseif(iflg.eq.2) then
        neq_test = n_awh_test
        neq_save = neq
        neq = 1
        n0 = neq
        i_awh1 = open_file('awh_mixed_variables.out','unknown')
        i_awh2 = open_file('awh_accum_terms.out','unknown')
        write(i_awh1,'(a)') 
     &  '** testing water,ngas,heat; ngas mixed fluid properties**'
        write(i_awh2,'(a)') 
     &  '** testing water,ngas,heat; mass energy accumulation terms**'
        do k =1, n_awh_test 
         phi_save = phi(k)
         phi(1) = phi_save
         pci_save = pci(k)
         pci(1) = pci_save
         t_save = t(k)
         t(1) = t_save
         s_save = s(k)
         s(1) = s_save
         call varchk_AWH(1)        
         call thrmwc(0)
c mixed fluid properties
         write(i_awh1,502) phi(1)  
         write(i_awh1,503) pci(1),
     &   t(1),s(1),(var_h2o(i),i=4,6),(var_awh_param(i),i =1,4),
     &   (var_awh_param(i),i =9,12),(var_awh_param(i),i =5,8)
c mixed fluid accumulation terms
         write(i_awh2,504) phi(1)  
         write(i_awh2,503) pci(1),
     &   t(1),s(1),(var_awh_param(i),i =17,28)
         phi(k) = phi_save            
         pci(k) = pci_save
         t(k)  = t_save
         s(k) = s_save
        enddo
        neq = neq_save
c stop after comparing parameters
502   format('water pressure',f12.4,/,t8,'PC',t22,'T',t41,'S',t54,
     & 'rol',t66,'drolp',t80,'drolt',t93,'rov',t106,'drovp',t118,'drovt'
     & ,t133,'drovpc',t147,'enw',t163,'dhlp',t178,'dhlt',t193,'dhlpc',
     & t209,'env',t222,'dhvp',t236,'dhvt',t252,'dhvpc')           
503   format(1p,19(1x,g11.4)) 
504   format('water pressure',f12.4,/,t8,'PC',t22,'T',t41,'S',t51
     & ,'den',t65,'damp',t75,'dame',t92,'dacp',t103,'dene',t119,'daep'
     & ,t134,'daeh',t149,'daepc',t164,'denc',t180,'dacp',t196,'dach'
     & ,t212,'dacpc')

      if(iptty.ne.0) write(iptty,*)
     &  'AWH parameter comparison completed: stopping'
      if(iout.ne.0) write(iout,*)
     &  'AWH parameter comparison completed: stopping'
       stop
      endif
c
      return
      end

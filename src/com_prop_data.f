      module com_prop_data
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
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

c 
c gaz 052521 initial coding
c
      real*8, allocatable  ::  den_h2o(:,:)
      real*8, allocatable  ::  enth_h2o(:,:)
      real*8, allocatable  ::  visc_h2o(:,:)
      real*8, allocatable  ::  psat_h2o(:,:)
      real*8, allocatable  ::  humid_h2o(:,:)

      real*8, allocatable  ::  den_h2o_old(:,:)					   
      real*8, allocatable  ::  enth_h2o_old(:,:)				   
      real*8, allocatable  ::  visc_h2o_old(:,:)	
      real*8, allocatable  ::  xv_h2o(:)
          
      real*8, allocatable  ::  den_ngas(:,:)
      real*8, allocatable  ::  enth_ngas(:,:)
      real*8, allocatable  ::  visc_ngas(:,:)
      real*8, allocatable  ::  xnl_ngas(:,:)

      real*8, allocatable  ::  den_ngas_old(:,:)
      real*8, allocatable  ::  enth_ngas_old(:,:)
      real*8, allocatable  ::  visc_ngas_old(:,:)
      real*8, allocatable  ::  xnl_ngas_old(:,:) 
      
      real*8, allocatable  ::  phi_old(:)					   
      real*8, allocatable  ::  t_old(:)
      real*8, allocatable  ::  pci_old(:)
      
      integer, allocatable  ::  ieval_flag(:)	
      integer, allocatable  ::  ieos_old(:)	

      integer  itot_calls, ic_eval_ngas, itot_calls_ngas, ic_eval
      integer ieosd_old, ieosd_last
      real*8 pl_old, pcl_old, tl_old, pl_last, tl_last, pcl_last
      real*8 rol_last,drolp_last,drolt_last,rov_last,drovp_last
      real*8 enl_last,dhlp_last,dhlt_last,env_last,dhvp_last,dhvt_last
      real*8 xvisl_last,dvislp_last,dvislt_last,xvisv_last,dvisvp_last
      real*8 drovt_last, dvisvt_last, dxnlpc_last
      real*8 ratio_eval_tot
c gaz 121923 
      real*8 xnl_max, xnl_ini, xnl_chng_low, xnl_chng_high
      
c gaz 070821 moved p_tol,t_tol,pc_tol  to comai and control to nr_stop_ctr    

      real*8 roc,drocpc,droct,roc_last,drocpc_last,droct_last
      real*8 hcg,dhcgp,dhcgt,hcg_last,dhcgp_last,dhcgt_last 
      real*8 xvisc,dxviscp,dxvisct,xvisc_last,dxviscp_last,dxvisct_last
      real*8 xnl,dxnlp,dxnlt,xnl_last,dxnlp_last,dxnlt_last
c gaz 072821
      real*8 hcg_last4,hcg4,dhcgp_last5,dhcgp5,dhcgt_last6,dhcgt6   
      real*8 dxnlpc,dxviscpc
c gaz 112223 ctest made global (help in debugging)
      real*8 ctest 
      integer ihenryiso     
c gaz 071121    
      logical var_last_flag, var_old_flag
c gaz 093021
c derived from mao 2013 (programmed by M Osullivan)
c xnlm,value at: (20,0.1),(120,0.1),(170,0.1),(20,0.2),(120,0.2),(170,0.2),(20,1.0),(120,1.0),(170,1.0)  
c dxnlmt,value at: (20-120,0.1),(120-170,0.1),(20-120,0.2),(120-170,0.2),(20-120,1.0),(120-170,1.0)
c dxnlmp,value at: (20,0.1-0.2),(120,0.1-0.2),(170,0.1-0.2),(20,0.2-1.0),(120,0.2-1.0),(170,0.2-1.0)      
      real*8, dimension(3) :: tlm 
      real*8, dimension(3) :: plm 
      real*8, dimension(9) :: xnlm  
      real*8, dimension(6) :: dxnlmt 
      data dxnlmt / 1.26314E-05,7.4277E-07, 2.50003E-05, 1.4659E-06, 
     &        0.000114951, 6.56042E-06 /
      real*8, dimension(6) :: dxnlmp
      data dxnlmp / 0.01677311,0.004404202, 0.004042636, 0.015527899, 
     &        0.004284101,  0.003965694 /
      data tlm / 20., 120., 170. /
c      data plm / 0.1, 0.2, 1.0 /
      data plm / 1.d-4, 0.2, 1.0 /
c      data xnlm / 1.71E-03, 4.43E-04, 4.06E-04, 3.38E-03, 8.84E-04,
      data xnlm / 1.d-6, 4.43E-04, 4.06E-04, 3.38E-03, 8.84E-04,      
     &        8.10E-04, 1.58E-02, 4.31E-03, 3.98E-03 / 
      integer, dimension(3,3) :: imaox 
      data imaox / 1, 2, 3, 4, 5, 6, 7, 8, 9 /
      integer, dimension(3,3) :: imaodxt 
      data imaodxt/ 1, 1, 2, 3, 3, 4, 5, 5, 6 /   
      integer, dimension(3,3) :: imaodxp 
      data imaodxp  / 1, 1, 2, 3, 3, 4, 5, 5, 6 /      
c imaox(3,3) - xnlm position (calculated from tlm and plm indices)
c imaodxt(3,3) - dxnlmt position 
c imaodxtp(3,3) - dxnlmp position 
      
      end module com_prop_data

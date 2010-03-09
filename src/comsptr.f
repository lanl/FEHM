      module comsptr
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1 
!D1 PURPOSE
!D1 
!D1 Include file containing passed parameters and pointers related to
!D1 the streamline particle tracking option.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 Initial implementation: 30-APR-1998, Programmer: G. Roemer
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/comsptr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:00   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 N/A
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 None
!D4 
!***********************************************************************

      logical :: random_flag = .false.
      logical :: sptrak = .false.
      logical :: unstruct_sptr = .false.
      real*8, allocatable :: ddx(:),ddy(:),ddz(:)
      real*8, allocatable :: ggg(:,:)
      real*8, allocatable :: corn(:,:)
      real*8, allocatable :: ddxv(:)
      real*8, allocatable :: vx1(:)
      real*8, allocatable :: axv(:)
      real*8, allocatable :: x1(:),x2(:),x3(:)
      real*8, allocatable :: vx1bv(:),vx2bv(:) 
      real*8, allocatable :: dtx(:),x_new(:)
      real*8, allocatable :: ddyv(:)
      real*8, allocatable :: vy1(:)
      real*8, allocatable :: ayv(:)
      real*8, allocatable :: y1(:),y2(:),y3(:)
      real*8, allocatable :: vy1bv(:),vy2bv(:) 
      real*8, allocatable :: dty(:),y_new(:)
      real*8, allocatable :: ddzv(:)
      real*8, allocatable :: vz1(:)
      real*8, allocatable :: azv(:)
      real*8, allocatable :: z1(:),z2(:),z3(:)
      real*8, allocatable :: vz1bv(:),vz2bv(:) 
      real*8, allocatable :: dtz(:),z_new(:)
      real*8, allocatable :: tt(:),dt(:),x61(:),axyzm(:),ttt(:),ttp1(:) 
      integer, allocatable :: ioutt(:),istop(:),ijkv(:),ic_new(:)
      integer, allocatable :: irray(:,:),ijkvs(:)
      integer, allocatable :: exit_node(:)
      real*8, allocatable :: xo(:),yo(:),zo(:),ttpo(:)
      integer, allocatable :: ijkvss(:)
      real*8, allocatable :: vx(:),vy(:),vz(:)
      real*8, allocatable :: omega_partial(:)
      real*8, allocatable :: sigma_partial(:)
c     changed from al1, al2, ... to dispersivity1, dispersivity2,
c     to avoid conflict with other variable names in the code
      real*8, allocatable :: dispersivity1(:)
      real*8, allocatable :: dispersivity2(:)
      real*8, allocatable :: dispersivity3(:)
      real*8, allocatable :: dispersivity4(:)
      real*8, allocatable :: dispersivity6(:)
      real*8, allocatable :: dispersivity7(:)
      real*8, allocatable :: dispersivity8(:)
      real*8, allocatable :: vratio(:)
      real*8, allocatable :: divdwt(:)
      real*8, allocatable :: edtmax(:)
      real*8, allocatable :: gradd(:,:)
      integer, allocatable :: izonebtc(:,:)
      integer, allocatable :: zbtc(:)
      integer, allocatable :: totlast(:,:)
      integer, allocatable :: totalpart(:)
c zvd 24-Jul-07
c Added totalpart_ret, totlast_ret for importance sampling
c Changed totalpart back to an integer
      real*8, allocatable :: totalpart_ret(:)
      real*8, allocatable :: totlast_ret(:)

      integer, allocatable :: tprpflag(:)
      integer, allocatable :: part_id(:,:)

c...skelkar  7/23/01 ........................
      integer, allocatable :: izoneplum(:,:)
      integer, allocatable :: plum(:)
      integer, allocatable :: plumpart(:)
      integer nplum
c...........................................

      real*8, allocatable :: divd_omr(:,:)
      real*8, allocatable :: dpordx_omr(:,:)

      real*8  divd(3), dpordx(3), dtensor(3,3)
      real*8  divd_weight

      integer nzbtc
      integer num_part
      integer nsize_tprp
! ZVD - 05/05 Added for increasing btc time steps
      real*8 dtmx, dtmn
! ZVD - 06/05 Added for limiting path output (iprto > 0)
      integer linelim
      real*8 distlim, tplim
      integer iprt,iprto
      integer iprtr
      integer ial
      real*8 courant_factor
c      real*8 al1,al2,al3,al4,al5
      integer itm,ist
      integer nx,ny,nz
      real*8 x10,y10,z10
      real*8 xdim,ydim,zdim
      real*8 tt1
      integer nsptrprops
      parameter(nsptrprops=20)
      integer write_prop(nsptrprops)
      real*8 sptr_prop(nsptrprops)
      integer itensor
      integer irevers
! ZVD - 05/05 Added for increasing btc time steps
      integer delta_part, part_frac, part_steps
      real*8 part_mult, time_btc

c...skelkar  11/13/02 ........................
c if freez_time is gt.0, then ptrac3 is called only at the end of the 
c flow calculations, and for a velocity frozen at that time,  
c particle tracks are calculated for freez_time days
      real*8 freez_time
c....................................
! ZVD - 08/10/05 Added for saving particle loactions and time for restarting sptr runs
c s kelkar sep 17 05 - sptrx_flag
! Use ttp1 instead of time_sptr
!      real*8, allocatable :: time_sptr(:)
      
      logical btc_flag, alt_btc, pod_flag, xyz_flag, xyz_flag2  
      logical ip_flag, sptr_flag, sptr_snode, sptrx_flag, trans_flag
      logical sptr_exists
c
      real*8, allocatable :: well_radius(:) 

! ZVD - 08/11/05 Added for tracking time particle reaches a btc zone
      real*8, allocatable :: ttbtc(:,:)
! ZVD - 04/06/09 Added for output location of particles that have been removed from system dure to bad location
      integer, allocatable :: lastnode(:)
! ZVD - 07/15/09 Add to mark wheter or not the starting particle location has been output
      logical, allocatable :: pstart_out(:)
! ZVD - 02/16/2010 Add variable to flag whether particles outside the model domain should be included in the simualtion
      logical :: exclude_particle
! ZVD - 03/04/2010 Add variable to force output of particle location in sptr2 if its start time is greater than the simulation end time
      logical :: output_end

      end module comsptr

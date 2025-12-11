      subroutine check_rlp_carb
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
!D1
!D1 PURPOSE
!D1
!D1 Test calculation of relative permeabilities and capillary pressures
!D1 and derivatives and output values.
!D1
!***********************************************************************
      
      use comai, only : form_flag, idpdp, ierr, neq, nrlp, wdd
      use comrlp, only : ishisrlp, rlpnew, delta_sat, num_sat,    
     &     sat_out,rlp_group,nphases, rlp_phase
      use comci, only : rlf, rvf,drvef,drlef
      use comdti, only : n0
      use comco2
      use comdi, only : ieos, irlp, pcp, icap,ices,s
      use comcomp, only : ioil

      implicit none
      integer i, j, mi, ndummy, neqtemp,ido, k,jj
      integer, allocatable :: ieostemp(:)
      integer, allocatable :: irlptemp(:)
      integer, allocatable :: icaptemp(:)
      real*8, allocatable  :: fwtemp(:)
      real*8, allocatable  :: fltemp(:)
      real*8, allocatable  :: fgtemp(:)
      real*8, allocatable  :: pcptemp(:)
      real*8, allocatable  :: pcgtemp(:)
      real*8, allocatable  :: rlwtemp(:)
      real*8, allocatable  :: rlltemp(:)
      real*8, allocatable  :: rlvtemp(:)
      real*8, allocatable  :: icestemp(:)
      real*8, allocatable  :: drleftemp(:)
      real*8, allocatable  :: drveftemp(:) 
      real*8, allocatable  :: rvftemp(:)       
      real*8, allocatable  :: stemp(:)     
      real*8, allocatable  :: s2temp(:) 	         
      real*8 :: dum1, dum2, dum3
      character*100 form_string, title_string
      character labels(3)*15,pn*30
      data labels/'co2_liquid','co2_gas','3-phase'/
      j=nphases(2)
      neqtemp = neq
      if (idpdp .ne. 0) then
         allocate (fwtemp(2*neq), fltemp(2*neq), fgtemp(2*neq))
         allocate (ieostemp(2*neq), irlptemp(2*neq), icaptemp(2*neq))
         allocate (icestemp(2*neq))
         allocate (pcptemp(2*neq), pcgtemp(2*neq))
         allocate (rlwtemp(2*neq), rlltemp(2*neq), rlvtemp(2*neq))
      else
         allocate (fwtemp(neq), fltemp(neq), fgtemp(neq))
         allocate (ieostemp(neq), irlptemp(neq), icaptemp(neq))
         allocate (pcptemp(neq), pcgtemp(neq))
         allocate (rlwtemp(neq), rlltemp(neq), rlvtemp(neq))
	 allocate (icestemp(neq) , s2temp(neq))
	 allocate(drleftemp(neq), drveftemp(neq),stemp(neq),rvftemp(neq))
      end if
      ndummy = 0
      
      if(icarb>0) then
         fwtemp = fw
         fltemp = fl
         fgtemp = fg
         rlwtemp = rl_w
         rlltemp = rl_l
         rlvtemp = rl_v
      else
         stemp = s
         rvftemp = rvf
         drveftemp = drvef
         drleftemp = drlef
      endif
      
      ieostemp = ieos
      irlptemp = irlp
      icaptemp = icap
      icestemp=ices
      pcptemp = pcp
      pcgtemp = pcg

      if (num_sat .eq. 0) then
         neq = 1 / delta_sat + 1
      else
         neq = num_sat
      end if
      
      title_string = "Relative permeability and " //
     &     "Capillary pressure"
      if (form_flag .eq. 1) then
         form_string = 'variables = "sw" "sl" "sg" "rl_w" '//
     &        ' "rl_l" "rl_g" "cp"' 
         write(ishisrlp, '(a)') trim(form_string)
         write(ishisrlp, 230) 50., 95., trim(title_string)
         write(ishisrlp, 230) 50., 90., trim(wdd)
      else if (form_flag .eq. 2) then
         form_string = 'Saturation, Liquid, CO2, ' //
     &        'Capillary pressure'
         write(ishisrlp, '(a)') trim(form_string)
      else
         form_string = ' "sw" "sl" "sg" "rl_w" '//
     &        ' "rl_l" "rl_g" "cp"' 
         write(ishisrlp, '(a)') trim(title_string)
         write(ishisrlp, '(a)') trim(form_string)
      end if

c     store a range of saturations in fw or s
c     neq=neq-2
      do i = 1, neq		
c     
c     first handle single (uncoupled) phases
c     calculate wetting phase saturation
c     
         if(icarb>0.or.ioil>0) then
            if (num_sat .eq. 0) then
               fw(i) = (i - 1) * delta_sat
               fl(i) = 1. - fw(i)
               fg(i) = 1. - fw(i)
            else
               fw(i) = sat_out(i)
               fl(i) = 1. - fw(i)
               fg(i) = 1. - fw(i)
            end if
         else
            if (num_sat .eq. 0) then
               s(i) = (i - 1) * delta_sat
            else
               s(i) = sat_out(i)
            end if       		
         endif		   

      end do   
      s2temp=fw    
c     
      do j = 1, nrlp
         
         
c     loop over model numbers
         do i=1,neq 
            irlp(i) = j
            icap(i) = j
            if (idpdp .ne. 0) then
               irlp(i+neq) = j
               icap(i+neq) = j
            end if
         end do 
         
c loop over twice, once assuming co2 is in liquid phase, then gas
        ido=1
 103    fl=0;fg=0
      do i = 1, neq
         if (num_sat .eq. 0) then
            fw(i) = (i - 1) * delta_sat
         else
            fw(i) = sat_out(i)
         end if
           if(ido.eq.1) then   ! co2 liquid
            fl(i) = 1 - fw(i)
            ices(i)=1
           else   ! co2 gas
            fg(i) = 1 - fw(i)
            ices(i)=3
           endif
         ieos(i) = 2
         if (idpdp .ne. 0) then
c double porosity
            fw(i + neq) = fw(i)
            fl(i + neq) = fl(i)
            ieos(i + neq) = 2
         end if
      end do
         if (form_flag .eq. 1) then
c tecplot
c write out the model number
            write (ishisrlp, 240) j,labels(ido)
         else if (form_flag .eq. 2) then
c surfer for cvs
         write (ishisrlp, 240) j,labels(ido)
         else
         end if
         if (rlpnew) then
            call rlp_cap(0)
            if (idpdp .ne. 0) call rlp_cap(neq)
         else
               do mi = 1, neq
c     calculate multi-phase relative perms.
                  call rlperm_co2(0,0,mi,rl_w(mi),
     &                 drl_ww(mi),drl_wg(mi),rl_l(mi),drl_lw(mi),
     &                 drl_lg(mi),rl_v(mi),drl_vw(mi),drl_vg(mi))
c     calculate multi-phase cap. pres.
                  call rlperm_co2(0,1,mi,pcp(mi),
     &                 dpcpw(mi),dpcpg(mi),pcg(mi),dpcgw(mi),dpcgg(mi),
     &                 dum1,dum2,dum3)
               end do
            endif
c     write out the results
c            write (ishisrlp, 240) rlp_group(j), pn(k)
            do i = 1, neq
c               if(k.ne.26) then
c                  if(icarb==0.and.ioil==0) then
c                     write (ishisrlp, '(7(g16.9, 1x))') s(i), rlf(i),
c     &                    rvf(i),  0., pcp(i)                      
c                  else if(k==22) then
c                     write (ishisrlp, '(7(g16.9, 1x))') fl(i), rl_l(i),
c     &                    rl_v(i), 0., pcp(i)             
c                  else          
                     write (ishisrlp, '(7(g16.9, 1x))') fw(i), fl(i),
     &                    fg(i),rl_w(i),rl_l(i),rl_v(i),pcp(i)	
c                  end if
c               else
c     3-phase
c!     cap pressure is too complicated to report for 3-phase, so let's write out 0's.  
c                  write (ishisrlp, '(7(g16.9, 1x))') fw(i), 
c     &                 rl_w(i),rl_l(i),rl_v(i),0.		
c               endif        	   
            end do
c double porosity
         if (idpdp .ne. 0) then
            if (form_flag .eq. 1) then
               write (ishisrlp, 245) j
            else if (form_flag .eq. 2) then
            else
            end if
            do i = neq + 1, 2*neq
               write (ishisrlp, '(7(g16.9, 1x))') fw(i), rl_w(i),
     &              rl_l(i), rl_v(i),pcp(i) ,44.
c     &              rl_l(i), pcp(i), drl_ww(i), drl_lw(i), drl_lg(i)
            end do
         end if
        ido=ido+1
        if(ido.le.2) goto 103
      end do

      fw = fwtemp
      fl = fltemp
      fg = fgtemp
      ices=icestemp
      ieos=ieostemp
      irlp=irlptemp
      icap=icaptemp
      
      pcp = pcptemp
      pcg = pcgtemp

      neq = neqtemp

      deallocate (fwtemp, fltemp, fgtemp, ieostemp, irlptemp, icaptemp, 
     &     pcptemp, pcgtemp, rlwtemp, rlltemp, rlvtemp)
      close (ishisrlp)

 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')
 240  format ('zone t = "model ', i4,1x,a30, '"')
 245  format ('zone t = "model ', i4, ' matrix"')
c		stop
      end subroutine check_rlp_carb

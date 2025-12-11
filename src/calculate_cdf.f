      subroutine calculate_cdf(model_flag,ith,divs_ith,ret,cdf)
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
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
!D1 Find appropriate interval in the cdf for importance sampling weighting.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 27-Oct-06, Programmer: H. Viswanathan
!D2
!D2 $Log:   /pvcs.config/fehm90/src/calculate_cdf.f_a  $
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.6 Streamline particle-tracking module 
!D3
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************

      use comai, only : ierr
      use comdti
      use comsk
      implicit none
      integer i, model_flag, ith, divs_ith
      real*8 cdf, m, b, ret,  k_for, cdflog

      if(model_flag.eq.11.or.model_flag.eq.12) then
c     Hari find the appropriate interval in the cdf
         i=1
         do i=1,nprobdivs(1)
            if(ret.le.rcdiv(i,1)) goto 11111
         enddo
11111    if(i.gt.nprobdivs(1)) i=nprobdivs(1)
         if(i.gt.1) then
            m = (probdiv(i,1)-probdiv(i-1,1))/(rcdiv(i,1)-rcdiv(i-1,1))
            b= probdiv(i,1)-m*rcdiv(i,1)
            cdf = m*ret+b
         else
            cdf=probdiv(i,1)
         endif
         if(model_flag.eq.12) cdf = sqrt(cdf)

      elseif(model_flag.ge.13.and.model_flag.le.14) then
         if(ret .eq. 1.) then
            cdf = 0.0
         else if(ret .lt. 1.) then
            write(ierr,*) 'Stop in calculate_cdf'
            write(ierr,*) 'ret .lt. 1.  stop'
            stop
         else 
            k_for=(ret-1)*k_rev(divs_ith)
            cdflog=cint_kf(divs_ith)+slope_kf(divs_ith)*dlog10(k_for)
            cdf=10.**cdflog
            if(model_flag.eq.14.or.model_flag.eq.16) cdf = sqrt(cdf)
         end if
      endif

      return
      
      end subroutine calculate_cdf

c.....................................................................

      

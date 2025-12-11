      subroutine calculate_ret(model_flag,ith,divs_ith,cdf,ret)
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

c s kelkar April 11 07
      use comai, only : ierr
      use comdti
      use comsk
      implicit none
      integer i, model_flag, ith, divs_ith
      real*8 cdf, m, b, ret,  k_for, k_f_log, dprob

      if(model_flag.eq.12.or.model_flag.eq.14) then
         write(ierr,*)'cant use model_flag=12 or 14 with'
         write(ierr,*)'flag_method=1. stop in calculate_r'
         write(ierr,*)' called from impsample_ptrk. STOP'
         stop
      elseif(model_flag.eq.11) then
c     find the appropriate interval in the probdiv
         i=1
         do i=1,nprobdivs(1)
            if(cdf.le.probdiv(i,1)) goto 11111
         enddo
11111    if(i.gt.nprobdivs(1)) i=nprobdivs(1)
         if(i.gt.1) then
            dprob =(probdiv(i,1)-probdiv(i-1,1))
            if(dprob.gt.0.) then
               m =(rcdiv(i,1)-rcdiv(i-1,1))/dprob
               ret=rcdiv(i-1,1)+m*(cdf-probdiv(i-1,1))
            else
               write(ierr,*)' WARNING in calculate_ret.'
               write(ierr,*)' table cdf constant or decreasing.'
               write(ierr,*)' set ret=rcdiv(i-1,1) & CONTINUE'
               ret=rcdiv(i,1)
            endif
         else
            ret=rcdiv(i,1)
         endif
      elseif(model_flag.eq.13) then
         if(cdf .eq. 0.) then
            ret = 1.0
         else if(cdf .lt. 0.) then
            write(ierr,*) 'Stop in calculate_r'
            write(ierr,*) 'cdf .lt. 0.  stop'
            stop
         else 
            k_f_log=(dlog10(cdf)-cint_kf(divs_ith))/slope_kf(divs_ith)
            k_for=10.**k_f_log
            ret=1.+ k_for/k_rev(divs_ith)
         end if
      endif

      return
      
      end subroutine calculate_ret

c.....................................................................

      

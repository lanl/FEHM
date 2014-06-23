      subroutine generate_flow
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
!D1 Generate flow for the mixing model.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 17-SEP-04,  Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/generate_flow.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
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

      use comrtd, only : rtd_subdiv, time_subdiv, capf,
     2     cumi, nsubdiv, maxmix_flag, mean_residence_time, 
     3     mean_porosity, mean_fl_density
      use comflow, only : a_axy
      use comai, only : neq, neq_primary, gdpm_flag
      use combi, only : izonef, sx1, ka, igdpm, ngdpm_layers
      use comdi, only : sk, esk, pflow, wellim

      implicit none
      real*8 flow_rate
      integer imodel
      logical ini_file
      integer iprev
      character*4 macro
      integer one, a_axy_counter
      parameter(one=1)
      real*8 flowin
      real*8 flowinmin
      parameter(flowinmin = 1.d-20)
      real*8 sumvol
      integer i, j
      real*8 total_volume, previous_age, previous_capf, cum_volume
      real*8 fractional_volume, current_age, current_capf
      integer current_index
      real*8 tol_vol
      parameter (tol_vol=-1.d-10)

c     Compute a_axy terms, store in array.

c     First, determine constants used in developing flux array
      sumvol = 0.
      do i = 1, neq_primary
         sumvol = sumvol + sx1(i)
      end do
      flow_rate = sumvol*mean_porosity*mean_fl_density
     2     /mean_residence_time
c     Compute source/sink terms terms that would otherwise
c     have been read in from flow macro

c     for each grid cell, find current value of internal age
c     at the entrance to that cell by interpolation
c     then compute by difference the fraction of the total
c     flow that enters that cell from the side entrance
      ka(1) = -1
      wellim(1) = 1.e9
      pflow(1) = 1.
      total_volume = 0.
      do i = 1, neq_primary
         total_volume = total_volume + sx1(i)
      end do
      previous_age = 0.
      previous_capf = 0.
      cum_volume = 0.
      current_index = 1

      do i = 1, neq_primary
         fractional_volume = sx1(i)/total_volume
         cum_volume = cum_volume + fractional_volume
         find_age: do j = current_index, nsubdiv

            if(cumi(j)-cum_volume.gt.tol_vol.or.
     2           j.eq.nsubdiv) then

               if(j.eq.1) then
                  current_age = time_subdiv(j)*
     2                 (cumi(j)/cum_volume)
                  current_capf = capf(j)*
     2                 (cumi(j)/cum_volume)
                  flowin = maxmix_flag*flow_rate*
     2                 (current_capf-previous_capf)
               else
                  current_age = time_subdiv(j-1)+
     2                 (time_subdiv(j)-time_subdiv(j-1))*
     3                 (cum_volume-cumi(j-1))/(cumi(j)-cumi(j-1))
                  
                  current_capf = capf(j-1)+
     2                 (capf(j)-capf(j-1))*
     3                 (current_age-time_subdiv(j-1))/
     4                 (time_subdiv(j)-time_subdiv(j-1))

                  current_capf = min(1.d0,current_capf)
                  flowin = maxmix_flag*flow_rate*
     2                 (current_capf-previous_capf)
               end if

               if(abs(flowin).gt.flowinmin) then
                  esk(i) = -20.
                  sk(i) = flowin
                  ka(i) = 1
                  izonef(i) = 1
               end if
               current_index = j
               previous_age = current_age
               previous_capf = current_capf
               if(i.eq.1) then
                  a_axy(1) = -maxmix_flag*flow_rate
                  sk(1) = a_axy(1)
                  a_axy(2) = -a_axy(1)
                  a_axy_counter = 2
               elseif(i.eq.neq_primary) then
                  a_axy_counter = a_axy_counter + 1
                  a_axy(a_axy_counter) = -a_axy(a_axy_counter-1-iprev)
                  a_axy_counter = a_axy_counter + 1
                  a_axy(a_axy_counter) = sk(i)
               else
                  a_axy_counter = a_axy_counter + 1
                  a_axy(a_axy_counter) = -a_axy(a_axy_counter-1-iprev)
                  a_axy_counter = a_axy_counter + 1
                  a_axy(a_axy_counter) = sk(i)
                  a_axy_counter = a_axy_counter + 1
                  a_axy(a_axy_counter) = -(a_axy(a_axy_counter-2)+sk(i))
               end if
c     Allow for gdpm node connected to this primary node
               iprev = 0
               if(gdpm_flag.ne.0) then
                  imodel = igdpm(i)
                  if(ngdpm_layers(imodel).ne.0) then
                     a_axy_counter = a_axy_counter + 1
                     iprev = 1
                  end if
               end if

               exit find_age
            end if
         end do find_age
      end do

      end subroutine generate_flow

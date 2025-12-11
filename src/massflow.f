      subroutine massflow
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
C***********************************************************************
CD1 PURPOSE
CD1 
CD1 This subroutine calculates the node mass, flow, and retardination 
CD1 factors. 
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 09-JAN-96    S. Henderson   22      Add prolog.
CD2              B. Robinson            Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/massflow.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:20   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2 Instead of calculating the mass, flux, and retardation factor data in 
CD2 part_track, now they are precalculated to get rid of the if statements 
CD2 in part_track.   				APR 17, 1997  CLI
CD2
CD2    Rev 1.5   Wed Jan 10 10:54:44 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.4   Tue Jan 09 15:01:58 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   Mon Jan 08 10:43:42 1996   robinson
CD2 Algorithm no longer recomputes things once heat and mass is turned off
CD2 
CD2    Rev 1.2   08/07/95 11:48:54   awolf
CD2 Fixed for dpdp - loops indexing modified
CD2 
CD2    Rev 1.1   03/16/95 09:48:44   llt
CD2 added PVCS log history
CD2
CD2    Rev 1.0   03/16/95 09:00:44   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3  
CD3 2.3.5 Cell-based particle-tracking module
CD3  
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
CPS
CPS    Trak_type is 1 for tracking liquid, 2 for vapor. For either 
CPS    compute the mass and exiting flow rate for the current node.
CPS    f_m_box is the index of the a_axy or a_vxy flow rate array which
CPS    holds the term for fracture-matrix flow for dpdp problem
CPS
CPS       IF this is a liquid solute
CPS
CPS	    Begin of loop
CPS           Compute fluid mass in the cell
CPS
CPS           FOR each node connected to the current node
CPS             Add outlet mass flow rate to the total outlet flow rate
CPS             Store cumulative outlet mass flow rate value (new!!)
CPS           ENDFOR each node connected to the current node
CPS
CPS           Compute retardation factor for the cell
CPS           FOR each node connected to the current node
CPS             Normalize the cumulative mass flow rate by total. . .
CPS. . .        outlet flow rate
CPS	    End of loop
CPS
CPS       ELSEIF this is a vapor solute
CPS
CPS	    Begin of loop	
CPS           Compute fluid mass in the cell
CPS
CPS           FOR each node connected to the current node
CPS             Add outlet mass flow rate to the total outlet flow rate
CPS             Store cumulative outlet mass flow rate value (new!!)
CPS           ENDFOR each node connected to the current node
CPS
CPS           Compute retardation factor for the cell
CPS           
CPS           FOR each node connected to the current node
CPS             Normalize the cumulative mass flow rate by total. . .
CPS. . .        outlet flow rate
CPS         End of loop
CPS
CPS          Note - get rid of i3.gt.99 check.  Not needed
CPS
CPS       ENDIF
CPS
CPS************ END of block that should be made a subroutine
c*******************Start of cell_time************************ 
Chun added compart.h and removed real*8 p(*) 

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comflow
      use compart
      use davidi
      implicit none

      integer add_fact,f_m_box,i,i2,ith,jj,neq1, lbox,ntmp,ntmp1
      integer nsizep_sf, jj_frac

      real matpor
      real*8  flow_out, flow_in

      neq1=neq+1
      nsizep_sf=nelm(neq+1) - neq - 1
      if(.not.allocated(p_sf))then
         allocate(p_sf(nsizep_sf))
      end if

      if(iptty.ne.0) then
         write(iptty,*) 'Computing particle probability vector'
      end if
      if(trak_type(1).eq.1)then


         do i=1,neq
	    flow_in = 0.
            flow_out = 0.
            add_fact = 0
	    jj=itrc(i)
	    jj_frac=jj
	    lbox=i

          ntmp=nelm(lbox)-neq
          ntmp1=nelm(lbox+1)-neq1
          f_m_box=(nelm(neq1)-neq1)*2+lbox
          if (irdof .ne. 13 .or. ifree .ne. 0) then
             mass(i)=rolf(i)*ps(i)*s(i)*sx1(i)
          else
             mass(i)=rolf(i)*ps(i)*sx1(i)
          end if
          do i2=ntmp,ntmp1
            if(a_axy(i2).gt.0) then
              p(i2) = a_axy(i2)
              p_sf(i2) = a_axy(i2)
              flow_out = flow_out+a_axy(i2)
            else
              p(i2) = 0.
              p_sf(i2) = 0.
              flow_in = flow_in + a_axy(i2)
            end if
          enddo
c     Store total outflow from node, except for f/m flux
          flow_ottot(i) = flow_out


C                                      take care of flow between 
C                                      fracture and matrix in dpdp model
          if(idpdp.ne.0) then
             if(a_axy(f_m_box).gt.0) then
                flow_out=flow_out+a_axy(f_m_box)
             else
                flow_in=flow_in+a_axy(f_m_box)
             endif
          endif

c     Call routine to compute transport particle parameters

          call ptrk_trans_params

c     Normalize the p array

          call normalize_p

       end do


c     matrix nodes

       do i=neq+1,n0
	    flow_in = 0.
          flow_out = 0.
	    jj=itrc(i)

          lbox=i-neq
	    jj_frac=itrc(lbox)
          add_fact=nelm(neq+1)-neq1

          ntmp=nelm(lbox)-neq+add_fact 
          ntmp1=nelm(lbox+1)-neq1+add_fact 
          f_m_box=(nelm(neq1)-neq1)*2+lbox
          if (irdof .ne. 13 .or. ifree .ne. 0) then
             mass(i)=rolf(i)*ps(i)*s(i)*sx1(i)
          else
             mass(i)=rolf(i)*ps(i)*sx1(i)
          end if
          do i2=ntmp,ntmp1
            if(a_axy(i2).gt.0) then
              p(i2) = a_axy(i2)
              p_sf(i2-add_fact) = p_sf(i2-add_fact) + a_axy(i2)
              flow_out = flow_out+a_axy(i2)
            else
              p(i2) = 0.
c     would be putting p_sf contribution here also, but it's 0
              flow_in = flow_in + a_axy(i2)
            end if
          enddo
c     Get total flow out of the f/m pair
          flow_ottot(lbox) = flow_ottot(lbox)+flow_out

C                                      take care of flow between 
C                                      fracture and matrix in dpdp model
          if(idpdp.ne.0) then
             
             if(a_axy(f_m_box).gt.0) then
                flow_in=flow_in-a_axy(f_m_box)
             else
                flow_out=flow_out-a_axy(f_m_box)
             endif
          endif

c     Call routine to compute transport particle parameters

          call ptrk_trans_params

c     Normalize the p array

          call normalize_p

c     Normalize the p_sf array

          call normalize_p_sf

       end do

      else if (trak_type(1).eq.2) then
        do 200 i=1,n0
	  flow_in = 0.
          flow_out = 0.
          add_fact = 0
	  jj=itrc(i)
	  lbox=i
	  if(i.gt.neq)then
            lbox=i-neq
	    add_fact=nelm(neq+1)-neq1
	    jj_frac=itrc(lbox)
	  endif
          ntmp=nelm(lbox)-neq+add_fact 
          ntmp1=nelm(lbox+1)-neq1+add_fact 
          f_m_box=(nelm(neq1)-neq1)*2+lbox
          mass(i)=rovf(i)*ps(i)*(1.-s(i))*sx1(i) 
          do i2=ntmp,ntmp1
            if(a_vxy(i2).gt.0) then
              p(i2) = a_vxy(i2)
              flow_out = flow_out + a_vxy(i2)
            else
              p(i2) = 0.
              flow_in = flow_in + a_vxy(i2)
            end if
          enddo
          if(idpdp.ne.0) then
            if(add_fact.gt.0) then
              if(a_vxy(f_m_box).gt.0) then
                flow_in=flow_in-a_vxy(f_m_box)
              else
                flow_out=flow_out-a_vxy(f_m_box)
              endif
            else
              if(a_vxy(f_m_box).gt.0) then
                flow_out=flow_out+a_vxy(f_m_box)
              else
                flow_in=flow_in+a_vxy(f_m_box)
              endif
            endif
          endif

	  do ith=1,nspeci
	    if(diffflag(jj,ith).eq.0)then
	      matpor=max(1.d-30,ps(i))
	    else
	      matpor=max(1.e-30,matrix_por(jj))
	    endif
            Rf(i,ith)=denr(i)/matpor
            Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/((1.-
     +        s(i))*rovf(i))+1.
	  enddo
       
	  if( flow_out .ne. 0. ) then        
            if( flow_out .ge. -flow_in ) then
              flow_ot(i)= flow_out
            else
              flow_ot(i)= -flow_in
            end if
            p(ntmp) = p(ntmp)/flow_out
            do i2=ntmp+1,ntmp1
              p(i2)=p(i2-1)+ p(i2)/flow_out
            end do 
          else
            flow_ot(i)= 0.
            do i2=ntmp,ntmp1
              p(i2)=0.
            end do
          end if

200     continue
      endif


      return

      contains

************************************************************
      
      subroutine ptrk_trans_params
      implicit none
      do ith=1,nspeci
         if(diffflag(jj,ith).eq.0)then
            matpor=max(1.d-30,ps(i))
            Rf(i,ith)=denr(i)/matpor
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/(s(i)*rolf(i))+1.
            else
               Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/rolf(i)+1.
            end if
         else
            if(idpdp.eq.0) then
               matpor = matrix_por(jj)
               Rf(i,ith)=denr(i)/matpor
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/(s(i)*rolf(i))+1.
               else
                  Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/rolf(i)+1.
               end if
            else
               if(diffflag(jj_frac,ith).le.-5)then
                  if(i.le.neq)then
                     Rf(i,ith) =(rd_frac(jj,ith)+kcoll(jj,ith)*
     +                    rcoll(jj,ith))/(1.+kcoll(jj,ith))
                  else
                     matpor=max(1.d-30,ps(i))
                     Rf(i,ith)=denr(i)/matpor
                     if (irdof .ne. 13 .or. ifree .ne. 0) then
                        Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/(s(i)*rolf(i))+1.
                     else
                        Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/rolf(i)+1.
                     end if
                  end if
               else
                  if(i.le.neq) then
                     matpor = ps(i+neq)
                     matpor=max(1.e-30,matpor)
                     Rf(i,ith)=denr(i+neq)/matpor
                     if (irdof .ne. 13 .or. ifree .ne. 0) then
                        Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/(s(i+neq)*
     +                       rolf(i+neq))+1.
                     else
                        Rf(i,ith)=Rf(i,ith)*kd(jj,ith)/rolf(i+neq)+1.
                     end if
                else
                  if(iout.gt.0) then
                     write(iout,*) 'Fatal error'
                     write(iout,*) 'For a dpdp simulation,'
                     write(iout,*) 'Do not apply the matrix'
                     write(iout,*) 'diffusion particle tracking'
                     write(iout,*) 'to the matrix nodes, only'
                     write(iout,*) 'the fracture nodes'
                     stop
                  end if
                end if
             end if
          endif
       end if
      enddo
      return
      end subroutine ptrk_trans_params


************************************************************


      subroutine normalize_p
      implicit none

      if( flow_out .ne. 0. ) then        
         if( flow_out .ge. -flow_in ) then
            flow_ot(i)= flow_out
         else
            flow_ot(i)= -flow_in
         end if
         p(ntmp) = p(ntmp)/flow_out
         do i2=ntmp+1,ntmp1
            p(i2)=p(i2-1)+ p(i2)/flow_out
         end do 
      else
         flow_ot(i)= 0.
         do i2=ntmp,ntmp1
            p(i2)=0.
         end do
      end if

      
      return
      end subroutine normalize_p

************************************************************


      subroutine normalize_p_sf
      implicit none

      if( flow_ottot(lbox) .ne. 0. ) then        
         p_sf(ntmp-add_fact) = p_sf(ntmp-add_fact)/
     2        flow_ottot(lbox)
         do i2=ntmp-add_fact+1,ntmp1-add_fact
            p_sf(i2)=p_sf(i2-1)+ p_sf(i2)/flow_ottot(lbox)
         end do 
      else
         do i2=ntmp-add_fact,ntmp1-add_fact
            p_sf(i2)=0.
         end do
      end if

      
      return
      end subroutine normalize_p_sf

************************************************************

      end subroutine massflow

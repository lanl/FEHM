      subroutine cell_time( nelm, p, lbox, add_fact, flow_through, 
     1     rock_density,
     2     neq, matpor, porosity, liquid_density, gas_density,
     3     saturation,
     3     cell_volume, trak_type, kd_fluid, fluid_mass,
     4     retard_factor,idpdp)
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
CD1
CD1  PURPOSE
CD1
CD1  This routine uses the outlet mass flow rates of the given node
CD1  to calculate probabilities of a particle next entering a given
CD1  neighboring node.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/cell_time.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
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
CD3  2.3.5 Cell-based particle-tracking module
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
CPS************ Start of block that should be made a subroutine (cell_time) Note - not same as code, need to change code to reflect this!!!
CPS
CPS         IF this is a liquid solute
CPS
CPS           Compute fluid mass in the cell
CPS
CPS           FOR each node connected to the current node
CPS             Add outlet mass flow rate to the total outlet flow rate
CPS             Store cumulative outlet mass flow rate value (new!!)
CPS           ENDFOR each node connected to the current node
CPS
CPS           Compute retardation factor for the cell
CPS
CPS         ELSEIF this is a vapor solute
CPS
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
CPS           ENDFOR each node connected to the current node
CPS
CPS          Note - get rid of i3.gt.99 check.  Not needed
CPS
CPS         ENDIF
CPS
CPS************ END of block that should be made a subroutine
C***********************************************************************
      use comflow
      use comdti
      implicit none

      integer nelm(*),lbox,add_fact,idpdp,neq,i2,f_m_box
      real*8 p(*),flow_through,flow_out,flow_in,rock_density
      real*8 porosity,liquid_density,gas_density,saturation,cell_volume
      integer trak_type
      real kd_fluid,matpor
      real*8 fluid_mass,retard_factor

      retard_factor=rock_density/matpor
      flow_in = 0.
      flow_out = 0.
c     Trak_type is 1 for tracking liquid, 2 for vapor. For either 
c     compute the mass and exiting flow rate for the current node.
c     f_m_box is the index of the a_axy or a_vxy flow rate array which
c     holds the term for fracture-matrix flow for dpdp problem
      f_m_box=(nelm(neq+1)-neq-1)*2+lbox
      if (trak_type.eq.1) then
         fluid_mass=liquid_density*porosity*saturation*cell_volume
         do i2=(nelm(lbox)-neq)+add_fact,(nelm(lbox+1)-neq-1)+add_fact
            if (a_axy(i2).gt.0) then
               p(i2) = a_axy(i2)
               flow_out = flow_out+a_axy(i2)
            else
               p(i2) = 0.
               flow_in = flow_in + a_axy(i2)
            end if
         enddo
c     take care of flow between fracture and matrix for dpdp problems
         if (idpdp.ne.0) then
            if (add_fact.gt.0) then
               if (a_axy(f_m_box).gt.0) then
                  flow_in=flow_in-a_axy(f_m_box)
               else
                  flow_out=flow_out-a_axy(f_m_box)
               endif
            else
               if (a_axy(f_m_box).gt.0) then
                  flow_out=flow_out+a_axy(f_m_box)
               else
                  flow_in=flow_in+a_axy(f_m_box)
               endif
            endif
         endif
         retard_factor=retard_factor*kd_fluid/(saturation*
     +        liquid_density)+1.
         
      else if (trak_type.eq.2) then
         fluid_mass=gas_density*porosity*(1.-saturation)*
     +        cell_volume
         do i2=(nelm(lbox)-neq)+add_fact,(nelm(lbox+1)-neq-1)+add_fact
            if (a_vxy(i2).gt.0) then
               p(i2) = a_vxy(i2)
               flow_out = flow_out + a_vxy(i2)
            else
               p(i2) = 0.
               flow_in = flow_in + a_vxy(i2)
            end if
         enddo
         if (idpdp.ne.0) then
            if (add_fact.gt.0) then
               if (a_vxy(f_m_box).gt.0) then
                  flow_in=flow_in-a_vxy(f_m_box)
               else
                  flow_out=flow_out-a_vxy(f_m_box)
               endif
            else
               if (a_vxy(f_m_box).gt.0) then
                  flow_out=flow_out+a_vxy(f_m_box)
               else
                  flow_in=flow_in+a_vxy(f_m_box)
               endif
            endif
         endif
         retard_factor=retard_factor*kd_fluid/((1.-
     +        saturation)*gas_density)+1.
      endif
      if( flow_out .ne. 0. ) then
         
         if( flow_out .ge. -flow_in ) then
            flow_through = flow_out
         else
            flow_through = -flow_in
         end if
         
        p(nelm(lbox)-neq+add_fact) = p(nelm(lbox)-neq+add_fact)/flow_out
         do i2=(nelm(lbox)-neq+1)+add_fact,(nelm(lbox+1)-neq-1)+add_fact
            p(i2)=p(i2-1) + p(i2)/flow_out
         end do
         
      else
         flow_through = 0.
         do i2=(nelm(lbox)-neq)+add_fact,(nelm(lbox+1)-neq-1)+add_fact
            p(i2)=0.
         end do
      end if
      
      return
      end
      

      subroutine sx_combine_ani(iflg)
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
!D1  PURPOSE
!D1
!D1  This subroutine combines the dimensions of a variable, used in 
!D1  generating the finite element coefficients. It also breaks
!D1  connections when Boundary Conditions connect nodes.
!D1  This routine is for anisotropy                            
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/sx_combine.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  3.1     Finite Element Coefficient Generation
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!***********************************************************************

      use davidi
      use comii
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none
c
      real*8  xi,yi,zi,xkb,ykb,zkb,dis2
      real*8  area_dis 
      integer iflg,i,i1,i2,iz
      integer jj,kb,neqp1,icon
      integer icode,iface
      integer kb_adv_x,kb_adv_y,kb_adv_z
c
      neqp1=neq+1
      if(iflg.eq.1) then
c
c break connections for specified pressure nodes to other specified
c pressure nodes
c
c check for previously defined mdnodes
c don't do this for anisotropy            
c
c break connections if neccessary
c ka(i) lt 0 is pressure BC,ieos(i) lt 0  is large rservoir
c
c gaz 10-19-2001 can't all break connections with solid phase yet
c anisotropy requires we break connections on control volume boundaries
       if(ice.eq.0) then
        icon=0
        do i=1,neq
         if(ka(i).ne.0) then
          kb_adv_x = ncon_adv(i,1)
          kb_adv_y = ncon_adv(i,2)
          kb_adv_z = ncon_adv(i,3)
          if(kb_adv_x.ne.0)then
           if(ka(kb_adv_x).ne.0) then
           i1=icxani(i)+1 
           i2=icxani(i+1)
            do jj=i1,i2
              sx_x(jj-neqp1)=0.0d00
              icon=icon +1
            enddo
           endif
          endif
          if(kb_adv_y.ne.0)then
           if(ka(kb_adv_y).ne.0) then
           i1=icyani(i)+1 
           i2=icyani(i+1)
            do jj=i1,i2
              sx_y(jj-neqp1)=0.0d00
              icon=icon +1
            enddo
           endif 
          endif
          if(kb_adv_z.ne.0)then
           if(ka(kb_adv_z).ne.0) then
           i1=iczani(i)+1 
           i2=iczani(i+1)
            do jj=i1,i2
              sx_z(jj)=0.0d00
              icon=icon +1
            enddo
           endif
          endif
         endif
        enddo
c
c write out warning that connections have been broken
c
      if(icon.ne.0) then
       write(iout,*) 
     &  'Found BC to BC connection(s): ',icon,',(now set=0.0)'
       if(iptty.ne.0) then
       write(iptty,*) 
     &  'Found BC to BC connection(s): ',icon,',(now set=0.0)'
       endif
      endif
c
      endif
      endif
      return
      end

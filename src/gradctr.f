      subroutine gradctr(iflg)
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
CD1  This subroutine manages initization and boundary conditions 
CD1  that involve gradients  
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 19-APR-02, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/gradctr.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
**********************************************************************

      use comai
      use combi
      use comco2
      use comdi
      use comdti
      use comfi
      use commeth
      implicit none

      integer iflg,icode, izone, inode, idir, iroot,neq_total
      integer mi,neqp1,i,i1,i2,j,jj,ja,ngrad_max,i_2nd
      integer open_file,igradm, imodel, jk, node_second_0
c gaz 052019      
      integer, allocatable ::  izonef_tmp(:) 
      character*80 gradm_name
      character*5 dgradm
      character*80 gradmod_root
      logical grad_dum
      character*3 grad_all
      parameter(ngrad_max= 50)
      real*8 dist,var_inode 

c================================================================
      if(igrad.eq.0) return
c================================================================
      if(iflg.eq.0) then
c
c read input  
c cordg : reference coordinate
c idirg : gradient coordinate direction (1 = x,2 = y, 3 = z)    
c igradf=1 : initial pressure is a variable 
c igradf=2 : initial temperature is a variable 
c igradf=3 : initial saturation  is a variable 
c igradf=4 : specified pressure (or head) is a variable
c igradf=5 : specified inflow temperture (or enthalpy) is a variable
c xg1 is x coordinate at reference point
c yg1 is y coordinate at reference point
c zg1 is z coordinate at reference point
c var0 is variable at cordg
c grad1 is the linear gradient
c
         if(.not.allocated(izone_grad)) then
            allocate(izone_grad(ngrad_max))
            allocate(igradf(ngrad_max)) 
            allocate(cordg(ngrad_max))
            allocate(var0(ngrad_max))
            allocate(grad1(ngrad_max))
            allocate(idirg(ngrad_max))
            allocate(gradmod_filename(ngrad_max))
            allocate(igradmodnamlen(ngrad_max))
            allocate(igradmodelfile(ngrad_max))
            igradmd = 0
            gradmod_filename = ''
         end if
         iroot = 8
         gradmod_root(1:8)='gradtemp'
         read(inpt,*) ngrad
        do i = 1,ngrad
         grad_dum = .false.
         read(inpt,'(a80)') wdd1
          do jj = 1,78
           if(wdd1(jj:jj+2).eq.'all') grad_dum = .true.
          enddo
         backspace inpt
         if(grad_dum) then
          read(inpt,*) 
     &        grad_all,cordg(i),idirg(i),
     &        igradf(i),var0(i),grad1(i)
          izone_grad(i) = -1
         else
          read(inpt,*) 
     &        izone_grad(i),cordg(i),idirg(i),
     &        igradf(i),var0(i),grad1(i)
          if(izone_grad(i).ne.-1) then
           izone_grad(i) = abs(izone_grad(i))  
          endif 
         endif    
        enddo
c       
c     save zone file and parameters
c 
         gradm_name = ''
         igradmd = igradmd + 1   
         idgradmc = 10000+igradmd
         write(dgradm,'(i5)')idgradmc 
         gradm_name(1:iroot) = gradmod_root(1:iroot)
         gradm_name(iroot+1:iroot+1) ='.'
         gradm_name(iroot+2:iroot+5) = dgradm(2:5)
         gradm_name(iroot+6:iroot+11) = '.gradf'
         gradmod_filename(igradmd)(1:iroot+11) 
     &   = gradm_name(1:iroot+11)
         igradmodnamlen(igradmd)= iroot+11
c         
c complete name here
c
         igradmodelfile(igradmd) = open_file(gradm_name, 'unknown')    
         j = igradmodelfile(igradmd)
         write(j,'(a4, 1x, i9)') 'grad',ngrad       
         do i = 1, ngrad
            write(j,*) 
     &           izone_grad(i),cordg(i),idirg(i),
     &           igradf(i),var0(i),grad1(i)
         enddo
         write(j,'(a4)') 'end '
         write(j,'(8(1x,i9))') (izonef(i),i=1,n0)
         close (j)
      else if(iflg.eq.1) then 
c
c modify initial values and BC's
c
        do jj = 1,igradmd
c
c read gradient info from auxillary files 
c and zone list for that request
c    
         j = igradmodelfile(jj)
         i1=igradmodnamlen(jj)
         j = open_file (gradmod_filename(jj), 'old')
c         open(j,file=gradmod_filename(jj)(1:i1),
c     &    status ='unknown')
         
         read(j,*)dgradm(1:4) ,ngrad      
         do i = 1, ngrad
            read(j,*) 
     &           izone_grad(i),cordg(i),idirg(i),
     &           igradf(i),var0(i),grad1(i)
         enddo
c gaz 052019
         if(.not.allocated(izonef_tmp)) allocate(izonef_tmp(n0))
         read(j,'(a4)') dgradm(1:4)
         read(j,*) (izonef_tmp(i),i=1,n0)
c    
         if(iread.le.0) then
c code with no restart file is present     
          if(gdpm_flag.ge.3.and.gdpm_flag.le.6) then  
           neq_total = neq_primary
          else
           neq_total = neq
          endif          
            do izone=1,ngrad
               do inode=1,neq_total
                  if(izonef_tmp(inode).eq.izone_grad(izone)
     &                  .or.izone_grad(izone).eq.-1) then
                     idir = idirg(izone)
                     dist = cord(inode,idir)-cordg(izone)
                     if(igradf(izone).eq.1) then
                        pho(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.2) then
                        to(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.3) then
                        so(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.4) then
                        pflow(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.5) then
                        esk(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.-5) then
                        esk(inode)=-(var0(izone) + grad1(izone)*dist)
                     else if(igradf(izone).eq.6) then
                        phometh(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.7) then
                        pflowmeth(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.8) then
                        qflux(inode)=var0(izone) + grad1(izone)*dist
c RJP 04/10/07 added the following part for CO2
                     else if(igradf(izone).eq.9) then
                        phoco2(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.10) then
                        pflowco2(inode)=var0(izone) + grad1(izone)*dist
                     endif
                  endif
               enddo
c  calculate variables from gradients for gdpm and gdkm    
             if(gdpm_flag.ge.3.and.gdpm_flag.le.6) then
               do inode=1,neq_primary
                 imodel = igdpm(inode)
                 if(ngdpm_layers(imodel).ne.0) then
                  if(izonef_tmp(inode).eq.izone_grad(izone))then
c notice that the gdpm nodes are 1D   
                   node_second_0 = nelm(nelm(inode+1))             
                   do jk = 1, ngdpm_layers(imodel)
                     i_2nd = node_second_0 + (jk-1)
                     idir = idirg(izone)
                     dist =  cord(inode,1) - cord(i_2nd,1)
                     if(igradf(izone).eq.11) then
                        var_inode = pho(inode)
                        pho(i_2nd) = var_inode + grad1(izone)*dist
                     else if(igradf(izone).eq.12) then
                        var_inode = to(inode)
                        to(i_2nd)= var_inode + grad1(izone)*dist   
                     endif
                   enddo
                   endif
                  endif
               enddo                
             endif      
            enddo
         else
            do izone=1,ngrad
               do inode=1,n0
                  if(izonef_tmp(inode).eq.izone_grad(izone)) then
                     idir = idirg(izone)
                     dist = cord(inode,idir)-cordg(izone)
                     if(igradf(izone).eq.4) then
                        pflow(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.5) then
                        esk(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.-5) then
                        esk(inode)=-(var0(izone) + grad1(izone)*dist)
                     else if(igradf(izone).eq.7) then
                        pflowmeth(inode)=var0(izone) + grad1(izone)*dist
                     else if(igradf(izone).eq.8) then
                        qflux(inode)=var0(izone) + grad1(izone)*dist
c RJP 04/10/07 added following for CO2
                     else if(igradf(izone).eq.10) then
                        pflowco2(inode)=var0(izone) + grad1(izone)*dist
                     endif
                  endif
               enddo
            enddo
         endif
c        rewind j
c        write(j,*) 'gradctr file ', jj, ' read'
        close(j, status = 'delete')         
       enddo
         if(allocated(izone_grad)) then
            deallocate(izone_grad)
            deallocate(izonef_tmp)
            deallocate(igradf)
            deallocate(var0)
            deallocate(grad1)
            deallocate(cordg)
            deallocate(idirg)
            deallocate(gradmod_filename)
            deallocate(igradmodnamlen)
            deallocate(igradmodelfile)
         end if

      else if(iflg.eq.2) then
c


      else if(iflg.eq.3) then
c

      endif
c 
      return
      end                

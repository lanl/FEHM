      subroutine subsidence(iflg)
!***********************************************************************
! Copyright 2011. Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los
! Alamos National Laboratory (LANL), which is operated by Los Alamos
! National Security, LLC for the U.S. Department of Energy. The U.S.
! Government has rights to use, reproduce, and distribute this software.
! Neither the U.S. Government nor Los Alamos National Security, LLC or
! persons acting on their behalf, make any warranty, express or implied,
! or assumes any liability for the accuracy, completeness, or usefulness
! of the software, any information pertaining to the software, or
! represents that its use would not infringe privately owned rights.
!
! The software being licensed may be Export Controlled.   It may not be
! distributed or used by individuals or entities prohibited from having
! access to the software package, pursuant to United States export
! control laws and regulations.  An export control review and
! determination must be completed before LANS will provide access to the
! identified Software.
!***********************************************************************

c      
c apply freeze-thaw to porosity change  and gridblock thickness
c
      use comai
      use combi
      use comci, only : dil,dile,deqh,dqh,dqt,dq,enlf,delf,delef
      use comdi
      use comsi
      implicit none
      real*8 z_i,z_kb,delzi,delzkb,z_tol,dz_max,coeff_ps
      real*8 delx2,dely2,delz2,ratio,tl,cos_z2,dz_temp
      real*8 PorosCoeff,DeltaHeight,z_plot_i,DeltaPres,KozenyCoeff
      real*8 AvgPoros,Counter,SubsWaveLength,SubsTime
      real*8 KinViscInv,KinViscInv0,DKVIT,skold,qhold,permsd
      integer iflg,i,i1,i2,j,kb,jj,kk,kb_max,ToDrain
      integer open_file,UnitNum
      parameter (z_tol=1.d-30,KinViscInv0 = 1.d6)
      if(i_subsid.eq.0) return
      if(iflg.eq.0) then
c read input,perform initial calculations (identify bottom/top zones)  
       allocate(izone_bot(neq))
       allocate(izone_top(neq))
       allocate(izone_drain(neq))
       allocate(izone_main(neq))
       izone_bot = 0
       izone_top = 0
       izone_drain = 0
       izone_main = 0
       allocate(i_chk(neq))
       allocate(z_plot(neq))
       allocate(z_plot_old(neq))
       allocate(dzrg_sub(neq))
       z_plot_old = 0.0
       ! read in fixed (usually bottom) boundary zone(s)
       read(inpt,*) kk, (i_chk(i), i = 1,kk)
       do jj = 1, kk
        i1 = i_chk(jj)
        do i = 1,neq
         if(izonef(i).eq.i1) izone_bot(i) = 3
        enddo
       enddo
       ! read in top boundary zone(s)
       read(inpt,*) kk, (i_chk(i), i = 1,kk)
       do jj = 1, kk
        i1 = i_chk(jj)
        do i = 1,neq
         if(izonef(i).eq.i1) then
            izone_top(i) = 3
         end if
        enddo
       enddo
       ! read in drainage face
       read(inpt,*) kk, (i_chk(i), i = 1,kk)
       do jj = 1, kk
        i1 = i_chk(jj)
        do i = 1,neq
         if(izonef(i).eq.i1) then
            izone_drain(i) = 3
         end if
        enddo
       enddo
       ! read in ice-rich zone
       read(inpt,*) kk, (i_chk(i), i = 1,kk)
       do jj = 1, kk
        i1 = i_chk(jj)
        do i = 1,neq
         if(izonef(i).eq.i1) then 
            izone_main(i) = 3
         end if
        enddo
       enddo
      else if(iflg.eq.-1) then
c calculate dzrg if not already done
         call airctr(11,0)
         ! KCL allocate pnx0,pny0,pnz0 with neq
         ! set equal to pnx, pny, and pnz
         allocate(pnx0(neq))
         allocate(pny0(neq))
         allocate(pnz0(neq))
         allocate(MaxDrainPart(neq))
         allocate(InitPoros(neq))
         allocate(InitPerm(neq))
         allocate(InitPres(neq))
         DrainInitFlag = 0
      else if(iflg.eq.1) then
c calculate subsidence
       if(itsat.eq.11) then
c first identify frozen and unfrozen soils
        ieos_aux = 0
        do i = 1, neq
         if(t(i).lt.IceEndTemp) then
          ieos_aux(i) = -1
         endif
        enddo 
       endif 
c 
c calculate changes in gridblock height
c
       if ((days>0).and.(DrainInitFlag==0)) then
          DrainInitFlag = 1
          do i = 1, neq
             pnx0(i) = pnx(i)
             pny0(i) = pny(i)
             pnz0(i) = pnz(i)
             psini(i) = ps(i)

c next line is for test only 
c             psini(i) = ((0.6-0.3)/1.e5)*(phi(i)-1.e5)+0.3

             MaxDrainPart(i) = 0.0
             InitPoros(i) = psini(i)
             InitPerm(i)  = pnz0(i)
             InitPres(i)  = phini(i)
          end do
       end if

       dzrg_sub = dzrg
       ToDrain = 1
       do i = 1, neq
         tl = t(i)
         coeff_ps  = PorosCoeff(i,tl,phi(i),ToDrain)
         ps(i) = psini(i)*coeff_ps

c for analytic comparision only!
c         ps(i) = 0.6*coeff_ps

         dzrg_sub(i) = dzrg(i) + 
     &      DeltaHeight(i,dzrg(i),ps(i))
         pnx(i) = KozenyCoeff(i,ps(i))*pnx0(i)
         pny(i) = KozenyCoeff(i,ps(i))*pny0(i)
         pnz(i) = KozenyCoeff(i,ps(i))*pnz0(i)

c for analytic comparision only!
c         pnx(i) = 1.e-7 + 5.e-7*(phi(i)-0.1D0)
c         pny(i) = pnx(i)
c         pnz(i) = pnx(i)
       enddo 
       
c       dzrg_sub = dzrg 
c       dzrg_sub(2460:2466) = 10.
c       dzrg_sub(2512:2561) = 10.
c       
c initialize z_plot 
c
       do i = 1,neq
        z_plot(i) = cord(i,igrav)
        i_chk(i) = 0
       enddo
       do jj = 1, neq
c start a fixed z displacement       
        if(izone_bot(jj).ne.0) then      
c        
c loop until all nodes are touched    
c need to first check for block centered FDM numerics
c for now this assumes 
c   
         if(ivf.eq.-1) then
          z_plot(jj) =  cord(jj,igrav) - dzrg_sub(jj)/2.
         endif
         i = jj
         
100      continue     
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         z_plot_i = z_plot(i)
         z_i = cord(i,igrav)
         delzi = dzrg(i)
         dz_max = 0.
         kb_max = i
         do j = i1,i2
          kb = nelm(j)
          if(kb.ne.i) then
           delx2=(cord(kb,1)-cord(i,1))**2
	       dely2=(cord(kb,2)-cord(i,2))**2
	       if(icnl.eq.0) dely2 = dely2 + (cord(kb,3)-cord(i,3))**2	      
           z_kb = cord(kb,igrav)-z_tol
           cos_z2 = (z_kb-z_i)**2/(delx2+dely2)
           dz_temp = (z_kb-z_i)*cos_z2
           if(dz_temp.gt.dz_max) then
            dz_max = dz_temp
            kb_max = kb
           endif
          endif
         enddo
         if(kb_max.ne.i) then
c
c  calculate a new z coordinate
c
          ratio = (dzrg_sub(i)+dzrg_sub(kb_max))/(dzrg(i)+dzrg(kb_max))
          z_plot(kb_max) = z_plot_i + dz_max*ratio
          i = kb_max
          go to 100
         else
          go to 101
         endif
        endif
101    continue        
       enddo   
c check for completeness    
       do i = 1, neq
        if(izone_bot(i).ne.0) then
        endif
       enddo
c KCL inserted to convert from coords to displacements
       do i = 1, neq
          z_plot(i) = z_plot(i) - cord(i,igrav)
       end do
      elseif (iflg.eq.2) then
         UnitNum = open_file('Predict.out','unknown')

         ! calculate and write average porosity in icy layer
         Counter = 0.0
         AvgPoros = 0.0
         do i = 1,neq
            if ((izone_main(i)==3).or.(izone_top(i)==3)) then
               AvgPoros = AvgPoros + ps(i)
               Counter = Counter + 1.0
            end if
         end do
         AvgPoros = AvgPoros/Counter
         write(UnitNum,*) 'avg porosity: ',AvgPoros
         write(*,*) 'avg porosity: ',AvgPoros

         ! calculate and write subsidence wavelength
         call GetSubsWaveLength(SubsWaveLength)
         write(UnitNum,*) 'subsidence wavelength: ',SubsWaveLength
         write(*,*) 'subsidence wavelength: ', SubsWaveLength

         close(UnitNum)
      else if(iflg.eq.3) then
c
c modify source/sink term for ice
c
c        return
        i = iad
        if(l.eq.0) return
        do i = 1, neq
          if(izone_drain(i).ne.0) then
            KinViscInv =  dil(i)
            
c            ratio = KinViscInv/KinViscInv0
            DKVIT = dile(i)/KinViscInv0
            if(t(i).lt.IceEndTemp) then
             ratio = 0.0
             DKVIT = 0.0
            else if(t(i).gt.LiqEndTemp) then
             ratio = 1.0
             DKVIT = 0.0
            else
             ratio = ((t(i)-IceEndTemp)/(LiqEndTemp-IceEndTemp))**2
             DKVIT = 2.0*((t(i)-IceEndTemp)/(LiqEndTemp-IceEndTemp))/
     &        (LiqEndTemp-IceEndTemp)
            endif
            permsd=abs(wellim(i))
            sk(i) = permsd*ratio*(phi(i)-pflow(i))     
            dq(i) = permsd*ratio
            dqt(i) = permsd*DKVIT*(phi(i)-pflow(i))
            if(sk(i).ge.0.0) then
             qh(i) = sk(i)*enlf(i)
             dqh(i) = dq(i)*enlf(i)+sk(i)*delf(i)
             deqh(i) = dqt(i)*enlf(i) + sk(i)*delef(i)
            else
             qh(i) = sk(i)*eflow(i)
             dqh(i) = dq(i)*eflow(i)
             deqh(i) = dqt(i)*eflow(i)           
            endif
          endif
        enddo
             
      end if

      end

      subroutine set_mptr(itmpin,failed_nodes,ith,istp,ithp,
     2     idaughter) 
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
CD1  This subroutine sets up the starting nodes and times for particle 
CD1  tracking, activated with macro mptr.  It is called only once, the 
CD1  first time the particle tracker is invoked.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              S. Henderson           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/set_mptr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2 modified to accommodate the using of variable confactors in 
CD2 rip-fehm runs.                                  2/6/98.
CD2 Add changes to accommodate the ingrowth model.
CD2                                                 July 19, 1997 cli
CD2
CD2 Made changes for multiple sepcies simulations and transient source
CD2 terms.                Apr. 18, 1997 cli
CD2
CD2    Rev 1.8   Wed Mar 27 11:35:48 1996   hend
CD2 Updated for ptrk output options
CD2 
CD2    Rev 1.7   Tue Mar 05 16:20:20 1996   hend
CD2 Added initialization for p array
CD2 
CD2    Rev 1.6   Fri Jan 12 17:56:26 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.5   Thu Jan 11 10:25:56 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   08/07/95 11:39:40   awolf
CD2 Fixed for DPDP (n0 instead of neq). Also changed for 
CD2 breakthrough output
CD2 
CD2    Rev 1.3   05/08/95 12:49:10   robinson
CD2 Calls plotc1 to write coordinate info to .trc file
CD2 
CD2    Rev 1.2   03/15/95 17:05:34   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.1   02/02/95 15:22:50   llt
CD2 added pvcs log info
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
c 
c Date: Summer '94
c Author: Stephen Henderson, EES-5
c         Contact: Bruce Robinson, EES-4
c                  George Zyvoloski, EES-5
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comflow
      use compart
      use davidi
      implicit none

      integer:: i,i2,ithp,num,split,itmp,ith,istp,j,k,lns2
      integer:: iold,iold1,jj,lns,lns1,kin1,itmpin,kin2,kin2_1,rloc
      integer:: idaughter,lstart,nps, failed_nodes !cli_added failed_nodes
      real:: ran_sp
      real*8:: inc,inflow,tmim, denom

      real, parameter:: dayinsecs=86400. 

 !cli_added 8/11/1999

c Set the initial node position and starting times for every particle. The
c manner in which these are set is dependent on the type of problem being 
c run. Pcnsk must either be all positive or all negative. Positive 
c indicates injecting according to source flow, while negative indicates
c evenly distributing the particles over the nodes. Pcnsk gives the 
c relative factor to distribute -- i.e. a 2 indicates twice as many 
c particles in this node or set of nodes. 
c make sure pcnsk is either always positive or always negative
         
 
 !set up failed package nodes for radionuclide release. 
 !failed_nodes<0,for single node release; otherwise,group of nodes  
      iold=ithp
      if(failed_nodes < 0 )then
         kin1=-failed_nodes
         kin2=-failed_nodes 
         kin2_1=1
      else
         kin1=insnode(itmpin-1)
         kin2=kin1+failed_nodes
         kin1=kin1+1
         kin2_1=failed_nodes
      endif

 !inject particles into the system

      if(pcnsk(kin1) > 0.) then
   !constant concentration release
         inflow=0
         if((trak_type(ith).eq.1).and.(ico2.lt.0)) then
  !Liquid only and isothermal air-water
            do i=kin1,kin2
               j=ptindex(i)
               if(sk(j).lt.0.) then
                  inflow=inflow+abs(sk(j))*pcnsk(i)
               endif
            enddo
            do k=1,istp
               tmim=t2sk(k)
               itmp=-1
               do i=kin1,kin2
                  j=ptindex(i)
                  if(sk(j).lt.0.)then
                     num=int(dum_p(k)*abs(sk(j))*pcnsk(i)/inflow)
                     if(num.ne.0)then
                        inc=(t2sk(k)-t1sk(k))/num
                        itmp=1
                        if(inc.lt.tmim)tmim=inc
                        do i2=1,num
                           ithp=ithp+1
                           box(ithp,ith)=j
                           frac_done(ithp,ith)=0.
                           theta(ithp,ith)=-99999.
                           start_time(ithp,ith)=
     2                          dayinsecs*(t1sk(k)+inc*(i2-1))
                        enddo
                     endif
                  endif
               enddo
               if(itmp.eq.1)then
	            !converte 2-D arrays into 1-D arrays 1/26/05
                  nsegs(ith)=nsegs(ith)+1 
	            lns1=aidex(ith)+nsegs(ith)
                  lns=lns1-1
	            nsport(lns1)=ithp
	            lstart=nsport(lns)+1
	            lsport(lns)=lstart 
	            nps=nsport(lns1)-nsport(lns)
                  dum_p(k)=nps
                  if(idaughter.ne.0)then
                     lns2=lns+lns
                     nprevd(lns)=0
                     tmsport(lns2-1)=t1sk(k)*dayinsecs
                     tmsport(lns2)=t2sk(k)*dayinsecs
                     ivdt(lns)=nps/(t2sk(k)-t1sk(k))
                     call sortranl(nps,start_time(lstart,ith), 
     2                    box(lstart,ith),tmim,rseed)
                  endif
               endif
            enddo
c     else if vapor only and isothermal air-water
         else if((trak_type(ith).eq.2).and.(ico2.lt.0)) then
            do i=kin1,kin2
               j=ptindex(i) 
               if(qh(j).lt.0.) then
                  inflow=inflow+abs(qh(j))*pcnsk(i)
               endif
            enddo
            do k=1,istp
               tmim=t2sk(k)
               itmp=-1
               do i=kin1,kin2
                  j=ptindex(i)
                  if(qh(j).lt.0.)then
                     num=int(dum_p(k)*abs(qh(j))*pcnsk(i)/inflow)
                     if(num.ne.0)then
                        inc=(t2sk(k)-t1sk(k))/num
                        itmp=1
                        if(inc.lt.tmim)tmim=inc
                        do i2=1,num
                           ithp=ithp+1
                           box(ithp,ith)=j
                           frac_done(ithp,ith)=0.
                           theta(ithp,ith)=-99999.
                           start_time(ithp,ith)=dayinsecs*(t1sk(k)+
     2                          inc*(i2-1))
                        enddo
                     endif
                  endif
               enddo
               if(itmp.eq.1)then
	            !converte 2-D arrays into 1-D arrays 1/26/05
                  nsegs(ith)=nsegs(ith)+1 
	            lns1=aidex(ith)+nsegs(ith)
                  lns=lns1-1
	            nsport(lns1)=ithp
	            lstart=nsport(lns)+1
	            lsport(lns)=lstart 
                  nps=nsport(lns1)-nsport(lns)
	            dum_p(k)=nps
                  if(idaughter.ne.0)then
                     lns2=lns+lns
                     nprevd(lns)=0
                     tmsport(lns2-1)=t1sk(k)*dayinsecs
                     tmsport(lns2)=t2sk(k)*dayinsecs
                     ivdt(lns)=nps/(t2sk(k)-t1sk(k))
                     call sortranl(nps,start_time(lstart,ith), 
     2                    box(lstart,ith),tmim,rseed)
                  endif
               endif
            enddo
c     else if liquid only and not isothermal air-water
         else if ((trak_type(ith).eq.1).and.(ico2.ge.0)) then
            do i=kin1,kin2
               j=ptindex(i)
               if((sk(j)*s(j).lt.0.)) then
                  inflow=inflow+abs(sk(j))*pcnsk(i)
               endif
            enddo
            do k=i,istp
               tmim=t2sk(k)
               itmp=-1
               do i=kin1,kin2
                  j=ptindex(i)
                  if(sk(j)*s(j).lt.0)then
                     num=int(dum_p(k)*abs(sk(j))*s(j)*pcnsk(i)/inflow)
                     if(num.ne.0)then
                        inc=(t2sk(k)-t1sk(k))/num
                        itmp=1
                        if(inc.lt.tmim)tmim=inc
                        do i2=1,num
                           ithp=ithp+1
                           box(ithp,ith)=j
                           frac_done(ithp,ith)=0.
                           theta(ithp,ith)=-99999.
                           start_time(ithp,ith)=dayinsecs*(t1sk(k)+
     2                          inc*(i2-1))
                        enddo
                     endif
                  endif
               enddo
               if(itmp.eq.1)then
	            !converte 2-D arrays into 1-D arrays 1/26/05
                  nsegs(ith)=nsegs(ith)+1 
	            lns1=aidex(ith)+nsegs(ith)
                  lns=lns1-1
	            nsport(lns1)=ithp
	            lstart=nsport(lns)+1
	            lsport(lns)=lstart 
                  nps=nsport(lns1)-nsport(lns)
	            dum_p(k)=nps
                  if(idaughter.ne.0)then
                     lns2=lns+lns
                     nprevd(lns)=0
                     tmsport(lns2-1)=t1sk(k)*dayinsecs
                     tmsport(lns2)=t2sk(k)*dayinsecs
                     ivdt(lns)=nps/(t2sk(k)-t1sk(k))
                     call sortranl(nps,start_time(lstart,ith), 
     2                    box(lstart,ith),tmim,rseed)
                  endif
               endif
            enddo
c     else if vapor only and not isothermal air-water
         else if((trak_type(ith).eq.2).and.(ico2.ge.0)) then
            do i=kin1,kin2
               j=ptindex(i)
               if(sk(j)*(1-s(j)).lt.0.) then
                  inflow=inflow+abs(sk(j))*pcnsk(i)
               endif
            enddo
            do k=1,istp
               tmim=t2sk(k)
               itmp=-1
               do i=kin1,kin2
                  j=ptindex(i)
                  if(sk(j)*(1-s(j)).lt.0.)then
                     num=int(dum_p(k)*abs(sk(j))*(1-s(j))*pcnsk(i)/
     2                    inflow)
                     if(num.ne.0)then
                        inc=(t2sk(k)-t1sk(k))/num
                        itmp=1
                        if(inc.lt.tmim)tmim=inc
                        do i2=1,num
                           ithp=ithp+1
                           box(ithp,ith)=j
                           frac_done(ithp,ith)=0.
                           theta(ithp,ith)=-99999.
                           start_time(ithp,ith)=dayinsecs*(t1sk(k)+
     2                          inc*(i2-1))
                        enddo
                     endif
                  endif
               enddo
               if(itmp.eq.1)then
	            !converte 2-D arrays into 1-D arrays 1/26/05
                  nsegs(ith)=nsegs(ith)+1 
	            lns1=aidex(ith)+nsegs(ith)
                  lns=lns1-1
	            nsport(lns1)=ithp
	            lstart=nsport(lns)+1
	            lsport(lns)=lstart 
                  nps=nsport(lns1)-nsport(lns)
	            dum_p(k)=nps
                  if(idaughter.ne.0)then
                     lns2=lns+lns
                     nprevd(lns)=0
                     tmsport(lns2-1)=t1sk(k)*dayinsecs
                     tmsport(lns2)=t2sk(k)*dayinsecs
                     ivdt(lns)=nps/(t2sk(k)-t1sk(k))
                     call sortranl(nps,start_time(lstart,ith), 
     2                    box(lstart,ith),tmim,rseed)
                  endif
               endif 
            enddo
    
         endif    

c     Take care of negative pcnsk
      else
            
         split=0
         do i=kin1,kin2
            split=split-pcnsk(i)
         enddo
         do k=1,istp
            tmim=t2sk(k)
            itmp=-1
            if(kin2_1.gt.dum_p(k))then
               do i=1,dum_p(k)
                  rloc=int(kin2_1*ran_sp(rseed))+kin1
                  ithp=ithp+1

CHari for fracture-matrix stuff
c                  write(6,*)'in set_mptr'
c                  write(6,*)'species','bin','fracflow'
c	            write(6,*)ith, itmpin-1,fracflow(ith,itmpin-1)
                  if(ripfehm.ne.1)then
                     j=ptindex(rloc)
                  else
                     if(ran_sp(rseed).le. fracflow(ith, itmpin-1))then
                        j=ptindex(rloc)
c                       write(6,*)'put particle in fracture',j
                     else
                        j=ptindex(rloc)+neq
c	                write(6,*)'put particle in matrix',j
                    endif
                 endif
                 box(ithp,ith)=j
                 frac_done(ithp,ith)=0.
                 theta(ithp,ith)=-99999.
                 start_time(ithp,ith)=dayinsecs*((t2sk(k)-t1sk(k))*
     2                ran_sp(rseed)+t1sk(k))
              enddo
              itmp=1
           else
              do i=kin1,kin2
CHari for fracture-matrix stuff
                 if(ripfehm.ne.1)then
                    j=ptindex(i)
                 else
                    if(ran_sp(rseed).le. fracflow(ith, itmpin-1))then
                       j=ptindex(i)
c	               write(6,*)'put particle in fracture',j
                    else
                       j=ptindex(i)+neq
c	               write(6,*)'put particle in matrix',j
                    endif
                 endif
                 num=int(dum_p(k)*(-pcnsk(i))/split+0.5)
                 if(num.ne.0)then
                    inc=(t2sk(k)-t1sk(k))/num
                    itmp=1
                    if(inc.lt.tmim)tmim=inc
                    do i2=1,num
                       ithp=ithp+1
                       box(ithp,ith)=j
                       frac_done(ithp,ith)=0.
                       theta(ithp,ith)=-99999.
                       start_time(ithp,ith)=dayinsecs*(t1sk(k)+
     2                      inc*(i2-1))
                    enddo
                 endif
              enddo
           endif
            if(itmp.eq.1)then
	        !converte 2-D arrays into 1-D arrays 1/26/05
              nsegs(ith)=nsegs(ith)+1 
	        lns1=aidex(ith)+nsegs(ith)
              lns=lns1-1
	        nsport(lns1)=ithp
	        lstart=nsport(lns)+1
	        lsport(lns)=lstart 
              nps=nsport(lns1)-nsport(lns)
	        dum_p(k)=nps
              if(idaughter.ne.0)then
                lns2=lns+lns
                nprevd(lns)=0
                tmsport(lns2-1)=t1sk(k)*dayinsecs
                tmsport(lns2)=t2sk(k)*dayinsecs
                ivdt(lns)=nps/(t2sk(k)-t1sk(k))
                call sortranl(nps,start_time(lstart,ith), 
     2               box(lstart,ith),tmim,rseed)
              endif
            endif
         enddo
  
      endif
 
      iold1=iold+1 
      num_particles(ith)=ithp

      if(pout.eq.0) then
         if(trak_type(ith).eq.1) then
            do i=iold1,ithp
               j=box(i,ith)
               jj=j+npt(ith)
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  denom = (rolf(j)*ps(j)*s(j)*sx1(j))
               else
                  denom = (rolf(j)*ps(j)*sx1(j))
               end if
               if(denom .ne. 0.) then
                  an(jj)=an(jj)+1/denom
               else
                  an(jj) = 0.
               end if
               anv(j)=an(j)
            enddo
         else
            do i=iold1,ithp
               j=box(i,ith)
               jj=j+npt(ith)
               denom = (rovf(j)*ps(j)*(1-s(j))*sx1(j))
               if(denom .ne. 0.) then
                  an(jj)=an(jj)+1/denom
               else
                  an(jj) = 0.
               end if
               anv(jj)=an(jj)
            enddo
         endif
      endif

      
      return 
      end






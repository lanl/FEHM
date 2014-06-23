      subroutine coneq1_gdpm(ndummy, matnum, spec_num)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To compute the Jacobian and residual terms of the concentration
CD1 equation, except for multiply defined node connections.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/coneq1.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2  
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.17   Wed May 29 14:31:34 1996   hend
CD2 Added variable diffusion with water content
CD2 
CD2    Rev 1.16   Fri May 24 09:52:12 1996   hend
CD2 Updated trac for mdnodes
CD2 
CD2    Rev 1.15   Fri May 03 14:19:20 1996   hend
CD2 Updated for GAZ mdnodes changes
CD2 Vel and Disp. set to 0 for parent connections
CD2 
CD2    Rev 1.14   Mon Apr 29 10:19:14 1996   hend
CD2 Updated to reflect GAZ changes -- all 3 sx components used
CD2 
CD2    Rev 1.13   Thu Apr 25 13:32:20 1996   hend
CD2 Updated for use in long/trans dispersion option
CD2 
CD2    Rev 1.12   Mon Mar 25 11:01:36 1996   hend
CD2 Fixed Indexing for Vapor Phase Velocities
CD2 
CD2    Rev 1.11   Thu Mar 21 13:17:42 1996   hend
CD2 Fixed to use indexing independent of istrw
CD2 
CD2    Rev 1.10   Mon Mar 04 16:06:28 1996   hend
CD2 Removed uneccessary calculations from coneq1 and added trac input option
CD2 
CD2    Rev 1.9   Mon Jan 29 13:55:40 1996   hend
CD2 Updated Requirements Traceability
CD2
CD2    Rev 1.8   12/13/95 10:27:32   robinson
CD2 Streamlined calculation of dispersion and advection terms
CD2 
CD2    Rev 1.7   09/29/95 16:12:16   llt
CD2 added small number on dividing -- was dividing by zero
CD2 
CD2    Rev 1.6   09/11/95 17:29:22   awolf
CD2 Dispersion elipsoid for 3-d calculations added
CD2 
CD2    Rev 1.6   08/29/95 08:41:08   awolf
CD2 Added elipsoidal representation of dispersivity for 3-D problems
CD2 
CD2    Rev 1.5   08/16/95 16:23:24   robinson
CD2 Corrected dispersion calculation for 0 saturation case
CD2 
CD2    Rev 1.4   08/14/95 16:15:42   awolf
CD2 Dispersion now allows (and requires) x,y,z components of alpha
CD2 
CD2    Rev 1.3   08/07/95 10:59:18   awolf
CD2 Dispersion now calculated with internode velocities.
CD2 Uses same algorithm as in flxo based on a_axy terms
CD2 
CD2    Rev 1.2   01/28/95 14:20:02   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.1   03/18/94 16:15:26   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:26   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name      Type        Description
CD3 
CD3 i          I          Node number of the current concentration
CD3                       unknown being worked on
c changed 8/5/94  to incorporate coupling
CD3 matnum       int         The current submatrix number
CD3 spec_num     int         The current coupled species number
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4 Global Subprograms
CD4
CD4 Name            Type     Description
CD4 ss_trans         N/A     calculate source/sink for tracer
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 npn, sx1, pnx, pny, pnz, phi, pcp, s, displx, disply, displz, 
CD4 dispvx, dispvy, dispvz, anl, anv, danl, danv, idualp, neq, nelm,
CD4 nelmdg, it8, it9, it10, it11, icnl, perml, permv, sx, igrav, grav,
CD4 t1, t2, t3, t4, t5, t5v, t6, t7, t8, t9, cord, dnwgta, upwgta, bp,
CD4 nrhs_sol, a, nmat_sol, awc, ayc, denci, dencj, rc, akc, drc, ps, 
CD4 istrw,rolf, rovf
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
C**********************************************************************
CD5 
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 itr          int         Solute unknown number
CD5 sx1d         real*8      Volume associated with the current node
CD5 dili         real*8      Liquid density
CD5 divi         real*8      Vapor density
CD5 anli         real*8      Liquid concentration
CD5 anlkb        real*8      Liquid concentration, connected node
CD5 anvi         real*8      Vapor concentration
CD5 anvkb        real*8      Vapor concentration
CD5 danli        real*8      Derivative of liquid concentration
CD5 danlkb       real*8      Derivative of liquid concentration
CD5 danvi        real*8      Derivative of vapor concentration
CD5 danvkb       real*8      Derivative of vapor concentration
CD5 anlri        real*8      Density times liquid concentration
CD5 anvri        real*8      Density times vapor concentration
CD5 danlri       real*8      Density times derivative of liquid
CD5                          concentration
CD5 danvri       real*8      Density times derivative of vapor
CD5                          concentration
CD5 icd          int         Integer index parameter
CD5 ii1          int         Integer index parameter
CD5 ii2          int         Integer index parameter
CD5 idg          int         Integer index parameter
CD5 iq           int         Integer index parameter
CD5 jmi          int         Integer index parameter
CD5 jml          int         Integer index parameter
CD5 jmia         int         Integer index parameter
CD5 jm           int         Do-loop index parameter
CD5 neqp1        int         Number of equations plus 1
CD5 ij           int         Do-loop index parameter
CD5 ij1          int         Integer index parameter
CD5 ij2          int         Integer index parameter
CD5 iz           int         Integer index parameter
CD5 kb           int         Integer index parameter
CD5 neighc       int         Integer index parameter
CD5 iw           int         Integer index parameter
CD5 sx2c         real*8      Real parameter used in calculation
CD5 sx3c         real*8      Real parameter used in calculation
CD5 sx4d         real*8      Real parameter used in calculation
CD5 sx4h         real*8      Real parameter used in calculation
CD5 sxzc         real*8      Real parameter used in calculation
CD5 radi         real*8      Parameter used in calculation
CD5 radkb        real*8      Parameter used in calculation
CD5 fid          real*8      Parameter used in calculation
CD5 fid1         real*8      Parameter used in calculation
CD5 axyd         real*8      Parameter used in calculation
CD5 axy          real*8      Parameter used in calculation
CD5 axyf         real*8      Parameter used in calculation
CD5 vxy          real*8      Parameter used in calculation
CD5 vxyf         real*8      Parameter used in calculation
CD5 dlaei        real*8      Derivative term used in calculation
CD5 dlaekb       real*8      Derivative term used in calculation
CD5 dvaei        real*8      Derivative term used in calculation
CD5 dvaekb       real*8      Derivative term used in calculation
CD5 heatc        real*8      Derivative term used in calculation
CD5 vxyd         real*8      Derivative term used in calculation
CD5 storage_term real*8      Derivative term used in calculation
CD5 iau          int         Integer index parameter
CD5 ial          int         Integer index parameter
CD5 kz           int         Integer index parameter
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following functions are carried out in this routine:
CD6 
CD6   Compute parameters needed later in the routine;
CD6   
CD6   In a loop over each node connected to the current node, the code
CD6   computes node numbers needed later in the calculation.
CD6   
CD6   The code next decides whether the calculation is two-dimensional
CD6   or three-dimensional, and carries out advection term
CD6   calculations accordingly, looping over each node connected to
CD6   the current node.
CD6   
CD6   The code next decides whether the solute is liquid or
CD6   vapor phase, and formulates the Jacobian terms and residual
CD6   terms by looping through each node connected to the current one,
CD6   and successively adding the corresponding advective or
CD6   dispersive term.
CD6   
CD6   Next, the code adds the storage and source/sink (or chemical
CD6   reaction) terms to the Jacobian and residuals.
CD6   
CD6   Finally, the code catches the special case of porosity of 0, and
CD6   sets terms accordingly before returning to the calling routine.
CD6     
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.3.4 Solute-transport equations
CD9 
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN coneq1
CPS 
CPS Compute terms used later in the calculation
CPS 
CPS IF this node is a matrix node
CPS   Set parameter used later to indicate matrix node
CPS ELSE
CPS   Set parameter used later to indicate primary node
CPS ENDIF
CPS 
CPS FOR each node connected to the current node
CPS 
CPS   Compute node numbers used later
CPS   
CPS   FOR each node connected to this node
CPS     IF this node is the current node
CPS       Set parameter accordingly
CPS     ENDIF
CPS     
CPS   ENDFOR
CPS   
CPS ENDFOR
CPS 
CPS IF this is a 3-D simulation
CPS 
CPS   FOR each node connected to the current node
CPS     Compute average permeability and dispersion coefficient
CPS   ENDFOR
CPS   
CPS ELSE this is a 2-D simulation
CPS 
CPS   FOR each node connected to the current node
CPS     Compute average permeability and dispersion coefficient
CPS   ENDFOR
CPS   
CPS ENDIF
CPS 
CPS IF this is a liquid phase or Henry's Law solute
CPS 
CPS   FOR each node connected to the current node
CPS     Compute liquid advection term
CPS   ENDFOR
CPS   
CPS 
CPS   FOR each node connected to the current node
CPS     Compute upwinding term
CPS   ENDFOR
CPS   
CPS 
CPS   FOR each node connected to the current node
CPS     Formulate Jacobian and residual terms of the current equation
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute dispersion term, add contribution to Jacobian and...
CPS     ... residual for liquid
CPS   ENDFOR
CPS ENDIF
CPS IF this is a vapor phase or Henry's law solute
CPS 
CPS   FOR each node connected to the current node
CPS     Compute vopor advection term
CPS   ENDFOR
CPS   
CPS 
CPS   FOR each node connected to the current node
CPS     Compute upwinding term
CPS   ENDFOR
CPS   
CPS 
CPS   FOR each node connected to the current node
CPS     Formulate Jacobian and residual terms of the current equation
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute dispersion term, add contribution to Jacobian and...
CPS     ... residual for vapor
CPS   ENDFOR
CPS   
CPS ENDIF
c changed 8/4/94 to allow for automatic adjustment of the rate 
c constants  (move to new routine react_add.f
CPS IF the porosity of this node is 0
CPS   Set Jacobian term to keep equation set well-behaved
CPS ENDIF
CPS 
CPS 
CPS END coneq1
CPS 
C**********************************************************************
c generate tracer equations
c
      use comflow
      use comcouple
      use comrxni
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comxi
      implicit none


      integer i,itr,icd,ii1,ii2,idg,iq,jmi,jml,jmia,jm,neqp1,ij
      integer ij1,ij2,iz,kb,neighc,iau,ial,kz, insp, kbnsp
      integer matnum,spec_num,ndummy,nmatavw,nmat2avw, nmatadd
      integer open_file, iwvfile
      real*8 sx1d,dili,divi,anli,anlkb,anvi,danli,danlkb,danvi
      real*8 sx3c,sxzc,anlri,anvri,danlri,danvri,sx2c,radi
      real*8 radkb,fid,fid1,axyd,axy,axyf,storage_term,dlaei
      real*8 dlaekb,heatc,vxyd,vxy,vxyf,anvkb,danvkb,dvaei
      real*8 t7wv(5000), t8wv(5000), t9wv(5000)
      real*8 dvaekb,toldil
      real*8 dvafi, dvafkb, anvi_wv, anvkb_wv, danvi_wv, danvkb_wv
      parameter(toldil = 1.d-20)
      real*8 slx2tt,svx2tt,tempx,tempy,tempz,templength,dispxavw
      real*8 dispzavw,alphaavw,axyflow,vxyflow,dilfkb,divfkb
      real*8 dilfi,divfi,dispyavw
      real*8 vmag,thetav,cord1x,cord1y,cord2x,cord2y,cord2xp,cord2yp
      real*8 cord1z,cord2z,newx,newy,newz,newdiff,concadiff,satr,satrkb
      real*8 rc_ss,drc_ss, pnxdum, pnydum, pnzdum

      real*8 newdiffkb, dumi, dumkb, dum_bar

      integer strindex,endindex,sehindexl,sehindexv
      integer pntr
      character*120 fname, root
      character*7 fsuffix
      integer iroot
      logical found_name, null1

      save t7wv, t8wv, t9wv, fname, found_name

c     Now called from ss_trans
c      call userc(1)
      if (ndummy.eq.0) then
         nmat_sol(1)=0
         nrhs_sol(1)=0
         strindex=1
         endindex=neq
         sehindexl=1
         sehindexv=1
      else if (ndummy.eq.neq) then
c         nmat_sol(1)=nmatb(4)
c         nrhs_sol(1)=nrhs(2)
         strindex=neq+1
         endindex=n0
         sehindexl=(sehsize/2)+1
         sehindexv=(sehsize/2)+1
      endif

      call dvacalc

      do i=strindex,endindex
         if (dispsame .eq. 0) then
            insp = i + (nsp - 1) * neq
         else
            insp = i
         end if
         itr=i+npn
         sx1d=sx1(i)
         dili=dil(i)
         anli=anl(itr)
         danli=danl(itr)
         anlri=anli*rolf(i)
         danlri=danli*rolf(i)
         if (irdof .ne. 13) then
            divi=div(i)
            anvi=anv(itr)
            danvi=danv(itr)
            anvri=anvi*rovf(i)
            danvri=danvi*rovf(i)
         else
            divi=0.            
            anvi=0.
            danvi=0.
            anvri=0.
            danvri=0.
         end if
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            satr = s(i)
         else
            satr = 1.d0
         end if

         neqp1=neq+1
         nmatavw = nelm(neqp1)-neqp1
         nmat2avw = 2*nmatavw
         if(i.gt.neq.and.idualp.eq.0) then
            icd=neq
            nmatadd = nmatavw
         else
            icd=0
            nmatadd = 0
         endif
         iz=i-icd

         ii1=nelm(iz)+1
         ii2=nelm(iz+1)
         idg=nelmdg(iz)-ii1+1
         iq=0
         jmi=nelmdg(iz)
         jmia=jmi-neqp1
         do jm=jmi+1,ii2
            iq=iq+1
            it10(iq)=istrw(jm-neqp1)
            if (sx(it10(iq),isox)+sx(it10(iq),isoy)+sx(it10(iq),isoz)
     &           .eq.0.) then
               iq=iq-1
            else
               it8(iq)=nelm(jm)+icd
               it9(iq)=jm-ii1+1
               it11(iq)=jm-neqp1
               ij1=nelm(nelm(jm))+1
               ij2=nelmdg(nelm(jm))-1
               do ij=ij1,ij2
                  if(nelm(ij).eq.iz) it12(iq)=ij-neqp1
               enddo
            endif
         enddo
         
         if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2) then

c-------------------------------------- PHS 9/2/2004
c   moving call for node i and calc of dumi to top of loop
c----------------------------------------------------
            if (dispsame .eq. 0) then
               newdiff = concadiff(1,mflagl(nsp,itrc(insp)),
     &              diffmfl(nsp,itrc(insp)),ps(i),satr,phi(i),t(i))
            else
               newdiff = concadiff(1,mflagl(1,itrcdsp(insp)),
     &              sehdiff(itrcdsp(insp)),ps(i),satr,phi(i),t(i))
            endif
            dumi  = satr*newdiff*ps(i)

c     
c     liquid phase calculations
c     
c     Find mass flux, upwind direction, compute term used later (liquid)
c
            do jm = 1, iq
               kb = it8(jm)
               iau=it11(jm)
               neighc=it9(jm)
               if(filter_flag(spec_num).eq.0) then
                 t8(neighc) = a_axy(iau+nmatadd)
               else
                 t8(neighc) = ftn_factor(istrw_cold(iau))*
     2                       a_axy(iau+nmatadd)
               endif
               fid = 0.5
               if(t8(neighc).lt.0.0) then
                  fid=dnwgta
                  fid1 = 1.-fid
                  t7(neighc) = fid*dil(kb)
                  t9(neighc) = fid1*dili
               elseif(t8(neighc).gt.0.) then
                  fid=upwgta
                  fid1 = 1.-fid
                  t7(neighc) = fid*dil(kb)
                  t9(neighc) = fid1*dili
               else
                  t7(neighc) = toldil
                  t9(neighc) = toldil
               end if
            enddo

c     
c     Add advection terms (liquid) to residual and Jacobian
c     
c     dir$ ivdep
            do jm=1,iq
               kb=it8(jm)
               kbnsp = kb + (nsp - 1) * neq
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               dilfi=t9(neighc)
               dilfkb=t7(neighc)
               anlkb=anl(kb+npn)
               danlkb=danl(kb+npn)
               axyflow = t8(neighc)
               axyd = axyflow/(dilfkb+dilfi)
               axyf=(dilfkb*anlkb+dilfi*anli)
               axy=axyd*axyf
               dlaei=axyd*dilfi*danli
               dlaekb=axyd*dilfkb*danlkb
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num))+axy
               bp(kz+nrhs_sol(spec_num))=bp(kz+nrhs_sol(spec_num))-axy
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+dlaei
               a(ial+nmat_sol(matnum))=a(ial+nmat_sol(matnum))-dlaei
               a(iau+nmat_sol(matnum))=a(iau+nmat_sol(matnum))+dlaekb
               a(jml+nmat_sol(matnum))=a(jml+nmat_sol(matnum))-dlaekb
            enddo

c
c     Formulate dispersion term (liquid) - different models for 
c     2 and 3 d (3D is icnl=0)
c     
            if (dispsame.eq.0) then
               if(icnl.eq.0) then
                  vmag=sqrt(pnx(n+i)*pnx(n+i)+
     +                 pny(n+i)*pny(n+i)+pnz(n+i)*pnz(n+i))
                  if (vmag.eq.0.) then
                     vmag=1e-30
                     pnxdum = 1e-30
                     pnydum = 0.
                     pnzdum = 0.
                  else
                     pnxdum = pnx(n+i)
                     pnydum = pny(n+i)
                     pnzdum = pnz(n+i)
                  end if
                  do jm=1,iq
                     kb=it8(jm)
                     kbnsp = kb + (nsp - 1) * neq
                     if (irdof .ne. 13) then
                        satrkb = s(kb)
                     else
                        satrkb = 1.0
                     end if
c------- PHS ---------- 9/3/2004 -----------------------------
                     newdiffkb = concadiff(1,mflagl(nsp,itrc(kbnsp)),
     &                    diffmfl(nsp,itrc(kbnsp)),ps(kb),satrkb,
     &                    phi(kb),t(kb))
               
                     dumkb = satrkb*newdiffkb*ps(kb)
                     dum_bar = 2*dumi*dumkb/(dumi+dumkb)
           if (i.le.neq_primary.and.kb.gt.neq_primary) then
	       dum_bar = dumkb
	     else if (i.gt.neq_primary.and.kb.le.neq_primary) then
	       dum_bar = dumi
           endif
c------- PHS ---------- 9/3/2004 -----------------------------

                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
c                    sx3c=sx(iw,2)
c                    sxzc=sx(iw,3)
                     if (ldsp.eq.0) then
                        cord1x=cord(iz,1)
                        cord1y=cord(iz,2)
                        cord1z=cord(iz,3)
                        cord2x=cord(kz,1)
                        cord2y=cord(kz,2)
                        cord2z=cord(kz,3)
                        dispxavw = tclx(nsp,itrc(insp))
                        dispyavw = tcly(nsp,itrc(insp))
                        dispzavw = tclz(nsp,itrc(insp))
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        tempz = (cord1z-cord2z)**2
                        templength=tempx+tempy+tempz
                        alphaavw=sqrt((dispxavw*dispyavw*dispzavw*
     +                       templength)/
     2                       (dispyavw*dispzavw*tempx +
     3                       dispxavw*dispzavw*tempy +
     4                       dispxavw*dispyavw*tempz) )
c------- PHS ---------- 9/3/2004 -----------------------------

                        slx2tt = sx2c*(sehvell(sehindexl)*alphaavw +
     &                       dum_bar)
c--------------------------------------------------------------------
                        sehindexl=sehindexl+1
                     else
                        call rotate(cord(iz,1),cord(iz,2),cord(iz,3),
     +                       pnxdum,pnydum,
     +                       pnzdum,
     +                       cord(kz,1),cord(kz,2),cord(kz,3)
     +                       ,newx,newy,newz)
                        cord1x=0.
                        cord1y=0.
                        cord1z=0.
                        cord2x=newx
                        cord2y=newy
                        cord2z=newz
                        dispxavw=tcly(nsp,itrc(insp))
                        dispyavw=tcly(nsp,itrc(insp))
                        dispzavw=tclx(nsp,itrc(insp))
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        tempz = (cord1z-cord2z)**2
                        templength=tempx+tempy+tempz
                        alphaavw = sqrt((dispxavw*dispyavw*dispzavw*
     +                       templength)/
     2                       (dispyavw*dispzavw*tempx +
     3                       dispxavw*dispzavw*tempy +
     4                       dispxavw*dispyavw*tempz) )
c------- PHS ---------- 9/3/2004 -----------------------------
                        slx2tt = sx2c*((vmag*alphaavw)+dum_bar)
c--------------------------------------------------------------------
                     endif
                     t5(neighc)=slx2tt
                  enddo
               else
                  radi=cord(iz,3)
                  vmag=sqrt(pnx(n+i)*pnx(n+i)+pny(n+i)*pny(n+i))
                  if (vmag.eq.0.) vmag=1e-30
                  thetav=acos(pnx(n+i)/vmag)
                  do jm=1,iq
                     kb=it8(jm)
                     kbnsp = kb + (nsp - 1) * neq
                     if (irdof .ne. 13) then
                        satrkb = s(kb)
                     else
                        satrkb = 1.0
                     end if
c------- PHS ---------- 9/3/2004 -----------------------------
                     newdiffkb = concadiff(1,mflagl(nsp,itrc(kbnsp)),
     &                    diffmfl(nsp,itrc(kbnsp)),ps(kb),satrkb,
     &                    phi(kb),t(kb))
                     dumkb = satrkb*newdiffkb*ps(kb)
                     dum_bar = 2*dumi*dumkb/(dumi+dumkb)
           if (i.le.neq_primary.and.kb.gt.neq_primary) then
	       dum_bar = dumkb
	     else if (i.gt.neq_primary.and.kb.le.neq_primary) then
	       dum_bar = dumi
           endif
c------- PHS ---------- 9/3/2004 -----------------------------
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     radkb=0.5*(radi+cord(kz,3))
                     sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
c                    sx3c=radkb*sx(iw,2)
                     if (ldsp.eq.0) then
                        cord1x=cord(iz,1)
                        cord1y=cord(iz,2)
                        cord2x=cord(kz,1)
                        cord2y=cord(kz,2)
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        dispxavw = tclx(nsp,itrc(insp))
                        dispyavw = tcly(nsp,itrc(insp))
                        alphaavw = sqrt((dispxavw*dispyavw * 
     +                       (tempx+tempy))/
     &                       (dispxavw*tempy + dispyavw*tempx ))
c----------------- PHS ------- 9/3/2004 -----------------------------
                        slx2tt = sx2c*(sehvell(sehindexl)*alphaavw +
     &                       dum_bar)
c--------------------------------------------------------------------
                        sehindexl=sehindexl+1
                     else
                        cord1x=0.
                        cord1y=0.
                        cord2xp=cord(kz,1)-cord(iz,1)
                        cord2yp=cord(kz,2)-cord(iz,2)
                        cord2x=
     +                       cos(thetav)*cord2xp-sin(thetav)*cord2yp
                        cord2y=
     +                       sin(thetav)*cord2xp+cos(thetav)*cord2yp
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        dispxavw = tclx(nsp,itrc(insp))
                        dispyavw = tcly(nsp,itrc(insp))
                        alphaavw = sqrt((dispxavw*dispyavw * 
     +                       (tempx+tempy))/
     &                       (dispxavw*tempy + dispyavw*tempx ))
c----------------------------- 9/3/2004 -----------------------------
                        slx2tt = sx2c*(vmag*alphaavw+dum_bar)
c--------------------------------------------------------------------
                     endif
                     t5(neighc)=slx2tt
                  enddo
               end if
            endif

c     
c     add dispersion term for liquid - 2 and 3d are the same 
c     at this point
c     
c     dir$ ivdep
            do jm=1,iq
               kb=it8(jm)
               kbnsp = kb + (nsp - 1) * neq
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               iw=it10(jm)
               if (dispsame.eq.1) then
                 if(filter_flag(spec_num).eq.0) then
                  heatc=sehvell(sehindexl)
                 else
                  heatc=ftn_factor(istrw_cold(iau))*sehvell(sehindexl)
                 endif
                  sehindexl=sehindexl+1
               else
                 if(filter_flag(spec_num).eq.0) then
                  heatc=t5(neighc)
                 else
                  heatc=ftn_factor(istrw_cold(iau))*t5(neighc)
                 endif
               endif
               axy=heatc*(anl(kb+npn)*rolf(kb)-anlri)
               dlaei=-heatc*danlri
               dlaekb=heatc*rolf(kb)*danl(kb+npn)
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num))+axy
               bp(kz+nrhs_sol(spec_num))=bp(kz+nrhs_sol(spec_num))-axy
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+dlaei
               a(ial+nmat_sol(matnum))=a(ial+nmat_sol(matnum))-dlaei
               a(iau+nmat_sol(matnum))=a(iau+nmat_sol(matnum))+dlaekb
               a(jml+nmat_sol(matnum))=a(jml+nmat_sol(matnum))-dlaekb
            enddo
         endif

**** Henry law addition ****
         if (icns(nsp).eq.-1.or.abs(icns(nsp)).eq.2) then

c----------------- PHS ------- 9/3/2004 -----------------------------
            if (dispsame .eq. 0) then
               newdiff = concadiff(2,mflagv(nsp,itrc(insp)),
     &              diffmfv(nsp,itrc(insp)),ps(i),satr,phi(i),t(i))
            else
               newdiff = concadiff(2,mflagv(1,itrcdsp(insp)),
     &              sehdiffv(itrcdsp(insp)),ps(i),satr,phi(i),t(i))
            endif
            dumi  = (1-satr) *newdiff*ps(i)
c-----------------------------------------------------------------

c     
c     vapor phase calculations
c     
c     Find mass flux, upwind direction, compute term used later (vapor)
c
            do jm = 1, iq
               kb = it8(jm)
               iau=it11(jm)
               neighc=it9(jm)
               t8(neighc) = a_vxy(iau+nmatadd)
               fid = 0.5
               if(t8(neighc).lt.0.0) then
                  fid=dnwgta
                  fid1 = 1.-fid
                  t7(neighc) = fid*div(kb)
                  t9(neighc) = fid1*divi
               elseif(t8(neighc).gt.0.) then
                  fid=upwgta
                  fid1 = 1.-fid
                  t7(neighc) = fid*div(kb)
                  t9(neighc) = fid1*divi
               else
                  t7(neighc) = toldil
                  t9(neighc) = toldil
               end if
            enddo
c     
c     Add advection terms (vapor) to residual and Jacobian
c     
c     dir$ ivdep
            do jm=1,iq
               kb=it8(jm)
               kbnsp = kb + (nsp - 1) * neq
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               divfi=t9(neighc)
               divfkb=t7(neighc)
               anvkb=anv(kb+npn)
               danvkb=danv(kb+npn)
               vxyflow = t8(neighc)
               vxyd = vxyflow/(divfkb+divfi)
               vxyf=(divfkb*anvkb+divfi*anvi)
               vxy=vxyd*vxyf
               dvaei=vxyd*divfi*danvi
               dvaekb=vxyd*divfkb*danvkb
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num))+vxy
               bp(kz+nrhs_sol(spec_num))=bp(kz+nrhs_sol(spec_num))-vxy
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+dvaei
               a(ial+nmat_sol(matnum))=a(ial+nmat_sol(matnum))-dvaei
               a(iau+nmat_sol(matnum))=a(iau+nmat_sol(matnum))+dvaekb
               a(jml+nmat_sol(matnum))=a(jml+nmat_sol(matnum))-dvaekb
            enddo
c
c     Formulate dispersion term (vapor) - different models for 
c     2 and 3 d (3D is icnl=0)
c     
            if (dispsame.eq.0) then
               if(icnl.eq.0) then
                  vmag=sqrt(pnx(n*2+i)*pnx(n*2+i)+
     +                 pny(n*2+i)*pny(n*2+i)+pnz(n*2+i)*pnz(n*2+i))
                  if (vmag.eq.0.) then
                     vmag=1e-30
                     pnxdum = 1e-30
                     pnydum = 0.
                     pnzdum = 0.
                  else
                     pnxdum = pnx(n*2+i)
                     pnydum = pny(n*2+i)
                     pnzdum = pnz(n*2+i)
                  end if
                  do jm=1,iq
                     kb=it8(jm)
                     kbnsp = kb + (nsp - 1) * neq
                     satrkb = s(kb)
c----------------- PHS ------- 9/3/2004 -----------------------------
                     newdiffkb = concadiff(2,mflagv(nsp,itrc(kbnsp)),
     &                    diffmfv(nsp,itrc(kbnsp)),ps(kb),satrkb,
     &                    phi(kb),t(kb))
                     dumkb = (1-satrkb)*newdiffkb*ps(kb)
                     dum_bar = 2*dumi*dumkb/(dumi+dumkb)
	     if (i.le.neq_primary.and.kb.gt.neq_primary) then
	       dum_bar = dumkb
	     else if (i.gt.neq_primary.and.kb.le.neq_primary) then
	       dum_bar = dumi
           endif
c-----------------------------------------------------------------
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     sx2c=sx(iw,isox) + sx(iw,isoy) + sx(iw,isoz)
c                    sx3c=sx(iw,2)
c                    sxzc=sx(iw,3)
                     if (ldsp.eq.0) then
                        cord1x=cord(iz,1)
                        cord1y=cord(iz,2)
                        cord1z=cord(iz,3)
                        cord2x=cord(kz,1)
                        cord2y=cord(kz,2)
                        cord2z=cord(kz,3)
                        dispxavw = tcvx(nsp,itrc(insp))
                        dispyavw = tcvy(nsp,itrc(insp))
                        dispzavw = tcvz(nsp,itrc(insp))
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        tempz = (cord1z-cord2z)**2
                        templength=tempx+tempy+tempz
                        alphaavw = sqrt((dispxavw*dispyavw*dispzavw*
     +                       templength)/
     2                       (dispyavw*dispzavw*tempx +
     3                       dispxavw*dispzavw*tempy +
     4                       dispxavw*dispyavw*tempz) )
c----------------- PHS ------- 9/3/2004 -----------------------------   
                        svx2tt = sx2c*(sehvelv(sehindexv)*alphaavw+
     &                        dum_bar)
c--------------------------------------------------------------------
                     else
                        call rotate(cord(iz,1),cord(iz,2),cord(iz,3),
     +                       pnxdum,pnydum,
     +                       pnzdum,
     +                       cord(kz,1),cord(kz,2),cord(kz,3)
     +                       ,newx,newy,newz)
                        cord1x=0.
                        cord1y=0.
                        cord1z=0.
                        cord2x=newx
                        cord2y=newy
                        cord2z=newz
                        dispxavw=tcvy(nsp,itrc(insp))
                        dispyavw=tcvy(nsp,itrc(insp))
                        dispzavw=tcvx(nsp,itrc(insp))
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        tempz = (cord1z-cord2z)**2
                        templength=tempx+tempy+tempz
                        alphaavw = sqrt((dispxavw*dispyavw*dispzavw*
     +                       templength)/
     2                       (dispyavw*dispzavw*tempx +
     3                       dispxavw*dispzavw*tempy +
     4                       dispxavw*dispyavw*tempz) )
c----------------- PHS ------- 9/3/2004 -----------------------------
                        svx2tt = sx2c * (vmag*alphaavw+dum_bar)
c--------------------------------------------------------------------
                     endif
                     t5v(neighc)=svx2tt            
                  enddo
               else
                  radi=cord(iz,3)
                  vmag=sqrt(pnx(n*2+i)*pnx(n*2+i)+
     +                 pny(n*2+i)*pny(n*2+i))
                  if (vmag.eq.0.) vmag=1e-30
                  thetav=acos(pnx(n*2+i)/vmag)
                  do jm=1,iq
                     kb=it8(jm)
                     kbnsp = kb + (nsp - 1) * neq
                     satrkb = s(kb)
c----------------- PHS ------- 9/3/2004 -----------------------------
                     newdiffkb = concadiff(2,mflagv(nsp,itrc(kbnsp)),
     &                    diffmfv(nsp,itrc(kbnsp)),ps(kb),satrkb,
     &                    phi(kb),t(kb))
                     dumkb = (1-satrkb)*newdiffkb*ps(kb)
                     dum_bar = 2*dumi*dumkb/(dumi+dumkb)
           if (i.le.neq_primary.and.kb.gt.neq_primary) then
	       dum_bar = dumkb
	     else if (i.gt.neq_primary.and.kb.le.neq_primary) then
	       dum_bar = dumi
           endif
c-----------------------------------------------------------------
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     radkb=0.5*(radi+cord(kz,3))
                     sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
c                    sx3c=radkb*sx(iw,2)
                     if (ldsp.eq.0) then
                        cord1x=cord(iz,1)
                        cord1y=cord(iz,2)
                        cord2x=cord(kz,1)
                        cord2y=cord(kz,2)
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        dispxavw = tcvx(nsp,itrc(insp))
                        dispyavw = tcvy(nsp,itrc(insp))
                        alphaavw =sqrt((dispxavw*dispyavw * 
     +                       (tempx+tempy))/
     &                       (dispxavw*tempy + dispyavw*tempx ))
c----------------- PHS ------- 9/3/2004 -----------------------------
                        svx2tt = sx2c*(sehvelv(sehindexv)*alphaavw+
     &                        dum_bar)
c--------------------------------------------------------------------
                        sehindexv=sehindexv+1
                     else
                        cord1x=0.
                        cord1y=0.
                        cord2xp=cord(kz,1)-cord(iz,1)
                        cord2yp=cord(kz,2)-cord(iz,2)
                        cord2x=
     +                       cos(thetav)*cord2xp-sin(thetav)*cord2yp
                        cord2y=
     +                       sin(thetav)*cord2xp+cos(thetav)*cord2yp
                        tempx = (cord1x-cord2x)**2
                        tempy = (cord1y-cord2y)**2
                        dispxavw = tcvx(nsp,itrc(insp))
                        dispyavw = tcvy(nsp,itrc(insp))
                        alphaavw =sqrt((dispxavw*dispyavw * 
     +                       (tempx+tempy))/
     &                       (dispxavw*tempy + dispyavw*tempx ))
c----------------- PHS ------- 9/3/2004 -----------------------------      
                        svx2tt = sx2c * (vmag*alphaavw+dum_bar)
c--------------------------------------------------------------------
                     endif
                     t5v(neighc)=svx2tt
                  enddo
               end if
            endif
c     
c     add dispersion term for vapor - 2 and 3D are the same 
c     at this point
c     
c     dir$ ivdep
            do jm=1,iq
               kb=it8(jm)
               kbnsp = kb + (nsp - 1) * neq
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               iw=it10(jm)
               if (dispsame.eq.0) then
                  heatc=t5v(neighc)
               else
                  heatc=sehvelv(sehindexv)
                  sehindexv=sehindexv+1
               endif
               vxy=heatc*(anv(kb+npn)*rovf(kb)-anvri)
               dvaei=-heatc*danvri
               dvaekb=heatc*rovf(kb)*danv(kb+npn)
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num))+vxy
               bp(kz+nrhs_sol(spec_num))=bp(kz+nrhs_sol(spec_num))-vxy
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+dvaei
               a(ial+nmat_sol(matnum))=a(ial+nmat_sol(matnum))-dvaei
               a(iau+nmat_sol(matnum))=a(iau+nmat_sol(matnum))+dvaekb
               a(jml+nmat_sol(matnum))=a(jml+nmat_sol(matnum))-dvaekb
c---------------------------------Temp spot for dumping i-kb wvflow
c------------------   only for nsp=1
c
               if(tort.GT.666) then
                if((i.EQ.1).AND.(kb.EQ.2)) then
                   if (.not. found_name) then
                      if (null1(root_name)) then
                         if (nmfil(5) .ne. nmfily(3) .and. 
     &                        nmfil(5) .ne. ' ') then
                            call file_prefix(nmfil(5), iroot)
                            if (iroot .gt. 100) iroot = 100
                            root(1:iroot) = nmfil(5)(1:iroot)
                         else
                            if (nmfil(2)(1:1) .eq. ' ' ) then
                               write (ierr, *) 
     &                             'FILE ERROR: nmfil2 file: ', 
     &                              nmfil(2), ' unable to determine ',
     &                              'wvflux file prefix'
                               stop
                            else
                               call file_prefix(nmfil(2), iroot)
                               if (iroot .gt. 100) iroot = 100
                               root(1:iroot) = nmfil(2)(1:iroot)
                            end if
                         end if                   
                      else
                         iroot = len_trim (root_name)
                         if (iroot .gt. 100) iroot = 100
                         root(1:iroot) = root_name(1:iroot)
                      end if
                      fsuffix = '.wvflux' 
                      fname =  root(1:iroot) // fsuffix
                      found_name = .true.
                   end if
                   iwvfile = open_file(fname, 'unknown')
                   rewind iwvfile
                   write(iwvfile,*) 'cord(i,1) nodei nodekb wvflow(kg/s)
     &                  wvflow(kg/yr) wlflow(kg/s) wlflow(kg/yr)'
                end if
c----     --    - -  --  -    - -   -     -       -        -

                pntr = iau+nmatadd

                if(i.LT.kb) then
                write(iwvfile,*) cord(i,1), i, kb, a_wvxy(pntr),
     &           a_wvxy(pntr)*3.15e7,a_axy(pntr), a_axy(pntr)*3.15e7
                end if
               end if
c----------------------------------------------------------------
            enddo
c--------------------------------------------------------------------
c Add in Water Vapor Diffusion term  for Henrys 3 (Water Vapor Tracer)
c    PHS 3/24/2004
c Find mass flux, upwind direction, compute term used later Watervapor
c
c    Take out t8 t7 t9 terms  replace with t7wv t8wv t9wv
c--------------------------------------------------------------------
           if(a_henry(nsp).EQ.666.) then
            do jm = 1, iq
               kb = it8(jm)
               iau=it11(jm)
               neighc=it9(jm)
               t8wv(neighc) = a_wvxy(iau+nmatadd)
               fid = 0.5
               if(t8wv(neighc).lt.0.0) then
                  fid=dnwgta
                  fid1 = 1.-fid
                  t7wv(neighc) = fid*dva(kb)
                  t9wv(neighc) = fid1*dva(i)
               elseif(t8wv(neighc).gt.0.) then
                  fid=upwgta
                  fid1 = 1.-fid
                  t7wv(neighc) = fid*dva(kb)
                  t9wv(neighc) = fid1*dva(i)
               else
                  t7wv(neighc) = toldil
                  t9wv(neighc) = toldil
               end if
            enddo
c
c     Add Water Vapor Diffusion terms  to residual and Jacobian
c       also, make anv = (ptot/pv)*anv
c
c     dir$ ivdep
            anvi_wv = anvi/(1-cnvf(i))
            danvi_wv = danvi/(1-cnvf(i))
            do jm=1,iq
               kb=it8(jm)
               kbnsp = kb + (nsp - 1) * neq
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               dvafi=t9wv(neighc)
               dvafkb=t7wv(neighc)
               anvkb_wv=anv(kb+npn)/(1-cnvf(kb))
               danvkb_wv=danv(kb+npn)/(1-cnvf(kb))

               vxyflow = t8wv(neighc)
               vxyd = vxyflow/(dvafkb+dvafi)
               vxyf=(dvafkb*anvkb_wv+dvafi*anvi_wv)
               vxy=vxyd*vxyf
               dvaei=vxyd*dvafi*danvi_wv
               dvaekb=vxyd*dvafkb*danvkb_wv
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num))+vxy
               bp(kz+nrhs_sol(spec_num))=bp(kz+nrhs_sol(spec_num))-vxy
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+dvaei
               a(ial+nmat_sol(matnum))=a(ial+nmat_sol(matnum))-dvaei
               a(iau+nmat_sol(matnum))=a(iau+nmat_sol(matnum))+dvaekb
               a(jml+nmat_sol(matnum))=a(jml+nmat_sol(matnum))-dvaekb
            enddo
           end if
c
c            END OF NEW WATER VAPOR SPECIES  PHS 3/24/2004
c--------------------------------------------------------------------
c
         endif
         call ss_trans(i,rc_ss,drc_ss)
         storage_term = sx1d*(awc*denci(itr)+ayc*dencj(itr))
         bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num)) +
     2        storage_term+rc_ss
         a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))+sx1d
     2        *awc*akc(itr)+drc_ss
      
      enddo

c------------------------- Temp set-up for dumping wvflow
      if(tort.EQ.666) close(iwvfile)
c---------------------------------------

      return
      end
      
      






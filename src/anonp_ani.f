      subroutine anonp_ani
!***********************************************************************
!  Copyright, 2004,  The  Regents of the University of California.
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
!D1 Categorize elements and call routines to generate finite element
!D1 coefficients. For anisotropic permeability.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    E!D
!D2 Date         Programmer     Number  Comments
!D2
!D2 09-01-2002   G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/anonp_ani.f_a  $
!D2 
!***********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.2 Control Volume Finite-Element Coefficient Generation (Anisotropic)
!D3
!***********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
!D4    FEHM Application Version 2.20
!D4
!***********************************************************************

      use combi
      use comdi
      use comei
      use comdti
      use comai
      implicit none

      real*8 al, alen1, alen2, alen3, det, difvol, dm(4,4), sumu
      real*8 tol, tp, tp1, tp2, tp3, tp4, vlorth, volest, volu
      real*8 vp, xie, xnelu, yie, ynelu, zero_d, zie, znelu 
      real*8 x1, x2, x3, x12, x13, x14, x15, x23, x34, x37, x45, x46 
      real*8 x56, x58, x67, x78, y1, y2, y3, y12, y13, y14, y15, y23 
      real*8 y34, y37, y45, y46, y56, y58, y67, y78, z1, z2, z3, z12
      real*8 z13, z14, z15, z23, z34, z37, z45, z46, z56, z58, z67, z78 
      integer i, i1, i2, i3, i4, i5, i6, ib, icode, iconn, id
      integer idiff, ie, ied, ig1, ij, ik, in, in1, in2, incon, int
      integer iortho, ipiv, is1, is2, isx, iu, j, je, jj, k, kb, kc, kj
      integer k1, k2, kk, lei, mm, n1, ncont, ned, neibr, nel, nele
      integer nelu, neluc, nelucd, neqp1, neu, nfd, nfu, node, nodeu
      integer nj, nsc, nsl, nslu, nt, numgb(8)
      integer ii, ipivkb, neq_total, n0_save, ipmax, lcnt, lcnt2

c using storage for a matrix
      real*8, allocatable ::  vole(:)
      real*8, allocatable ::  dum(:)
      integer, allocatable :: ncon(:)
      integer, allocatable :: iorth(:)
      integer, allocatable :: nf(:)
      integer, allocatable :: idum(:)
      integer, allocatable :: idum1(:)
      integer, allocatable :: idum2(:)
      integer, allocatable :: nodeuf(:)
      integer, allocatable :: inoduf(:)
      integer, allocatable :: narr(:,:)
      integer, allocatable :: nsf(:)
      integer, allocatable :: iplace(:)

      integer, allocatable :: ncon_elem(:,:)
      integer, allocatable :: ncon_pos(:,:)
c
      integer, allocatable :: ianni(:)
      integer, allocatable :: iannip(:,:)
      integer, dimension(12)  :: ione  
      integer, dimension(12)  :: itwo 
      integer, dimension(12)  :: ithree
      integer, dimension(12)  :: ifour   
      integer, dimension(12)  :: ifive   
      integer, dimension(12)  :: isix    
      integer, dimension(12)  :: iseven  
      integer, dimension(12)  :: ieight  
      integer, dimension(12)  :: inine   
      integer, dimension(12)  :: iten    
      integer, dimension(24)  :: indx    
      integer, dimension(8,8) :: ipos
      integer, dimension(8,3) :: iface
      integer, dimension(8,3) :: kb_face_local
      integer, dimension(8)   :: kkx
      integer, dimension(8)   :: kky
      integer, dimension(8)   :: kkz
   

      integer idani,kb1,kb2,irowc,minkb,maxkb
      integer kki
      integer max_elem,icx,icy,icz,ic,ncon_space,nedge
      integer ncon_space_sx,kb_node,kb_nodex,kb_nodey,kb_nodez
      integer min1kb,max1kb,min2kb,max2kb,ifaced 
      integer minkbx,maxkbx,minkby,maxkby,minkbz,maxkbz         
      integer itest_iso,j1,j2,ipos1,ipos2
      real*8 perm_scale,dis,cosx,cosy,cosz
      real*8 area_j,vol_con,sx_fac,termx,termy,termz,a_tol
      real*8 theta1,theta2,theta3,deg2rad           
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33
      real*8 ak11,ak12,ak13,ak22,ak23,ak33
      real*8  ,dimension(12)  :: xcell     
      real*8  ,dimension(12)  :: ycell     
      real*8  ,dimension(12)  :: zcell     
      real*8  ,dimension(12)  :: cell_dir     
      real*8  ,dimension(12)  :: area      
      real*8  ,dimension(24,24) :: a_ani 
      real*8  ,dimension(24,24) :: a_ani_c
      real*8  ,dimension(24,24) :: a_anii 
      real*8  ,dimension(8,6) :: aklem 

      real*8  ,allocatable :: sx_temp_x(:,:,:)
      real*8  ,allocatable :: sx_temp_y(:,:,:)
      real*8  ,allocatable :: sx_temp_z(:,:,:)

      real*8   ,allocatable :: idumx(:)
      real*8   ,allocatable :: idumy(:)
      real*8   ,allocatable :: idumz(:)
      integer   ,allocatable :: nelm_temp(:,:)
      logical iani_chk

      parameter (a_tol=1.d-15, deg2rad=0.017453293d00,iani_chk=.false.)

c  gaz modified ione and itwo

c      data ione   / 1,2,4,1,5,6,8,5,1,2,3,4 /  
c      data itwo   / 2,3,3,4,6,7,7,8,5,6,7,8 / 
      data ione   / 1,2,3,4,5,6,7,8,1,2,3,4 /  
      data itwo   / 2,3,4,1,6,7,8,5,5,6,7,8 / 

      data ithree / 1,6,2,5,4,7,3,8,9,10,11,12 / 

      data ifour  / 1,4,8,5,1,2,6,5,1,2,3,4 /   
      data ifive  / 2,3,7,6,4,3,7,8,5,6,7,8 /  
      data isix   / 1,1,1,1,0,0,0,0,0,0,0,0 /  
      data iseven / 0,0,0,0,1,1,1,1,0,0,0,0 /  
      data ieight / 0,0,0,0,0,0,0,0,1,1,1,1 /  
      data inine  / 2,1,4,3,6,5,8,7,10,9,12,11/
      data iten   / 4,3,2,1,8,7,6,5,12,11,10,9/
      data cell_dir / 1.,+1.,+1.,1.,1.,1.,1.,1.,1.,1.,1.,1./
      data neluc / 0 /
      data ipos /14,15,18,17,5,6,9,8,
     &           13,14,17,16,4,5,8,7,
     &           10,11,14,13,1,2,5,4,
     &           11,12,15,14,2,3,6,5,
     &           23,24,27,26,14,15,18,17,
     &           22,23,26,25,13,14,17,16,
     &           19,20,23,22,10,11,14,13,
     &           20,21,24,23,11,12,15,14/
      data iface /1,-1,-2,2,4,-4,-3,3,5,6,-6,-5,8,7,-7,-8,
     &      -9,-10,-11,-12,9,10,11,12/
      data kb_face_local /2,0,0,3,6,0,0,7,4,3,0,0,8,7,0,0,
     &      0,0,0,0,1,2,3,4/
      data kkx /1,0,0,2,3,0,0,4/
      data kky /1,2,0,0,3,4,0,0/
      data kkz /0,0,0,0,1,2,3,4/  


c     triple vector product
      tp(x1, y1, z1 ,x2, y2, z2, x3, y3, z3) = 
     .     x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 *z2 - 
     .     x3 * y2 * z1 - x2 * y1 * z3 - x1 * y3 *z2
c     length of vector 
      al(x1, y1, z1) = sqrt(x1 * x1 + y1 * y1 + z1 * z1)
c     vector product
      vp(x1, y1, x2, y2) = x1 * y2 - y1 * x2
      
      zero_d = 0.0
c
c GAZ 02-08-01
c set neq to neq_primary here for gdpm, then set back later
c also set n0 to neq_primary here for gdpm, then set back later
c
      if(gdpm_flag.ne.0 .or. nriver.ne.0) then
         neq_total = neq
         neq = neq_primary
         n0_save= n0
         n0 = neq_primary
      endif
      
      tol = 1.e-10
      
c     determine if orthogonal elements should be fully connected
c     they should always be connected for anisotropy
      iconn = 1
      
      allocate(iplace(n0))

       iplace=0
      
c
c calculate # connections for each node
c
      do ie=1,nei
         do j=1,ns
            i=nelm((ie-1)*ns+j)
            if (i .ne. 0) then
               iplace(i)=iplace(i)+1
            endif
         enddo
      enddo

      ipmax=0
c estimate maximum connections
      if(icnl.eq.0) then
       lcnt=8
       lcnt2=27
       ldn2=lcnt2*neq
      else
       lcnt=4
       lcnt2=9
       ldn2=lcnt2*neq
      endif
c add index space for ncon
c note "-2" is based on ns*connected elements vs difference stencils
      ldn2=ldn2+neq+1
      do i=1,neq
         ipmax=max(ipmax,iplace(i))
         if(iplace(i).gt.lcnt) then
          ldn2=ldn2+iplace(i)*(ns-2)+1-lcnt2
         endif
      enddo
      
      nemx=ipmax
      allocate(vole(nei),iorth(nei),nf(n0),narr(nei,8),nsf(nei))
      iorth=0

c     find number of nodes in each element
      do ij = 1, nei
         nsc = 0
         do ied = 1, ns
            ik = nelm((ij - 1) *ns + ied)
            if(ik .ne. 0) nsc = nsc + 1
            numgb(ied) = ik
            narr(ij, ied) = ied
         end do
         nsf(ij) = nsc
         
c     sort elements into ascending order
         do j = 2, ns
            n1 = numgb(j)
            do i = j - 1, 1, -1
               if(numgb(i) .lt. n1) go to 17
               narr(ij, i + 1) = narr(ij, i)
               numgb(i + 1) = numgb(i)
            end do
            i = 0
 17         narr(ij, i+1) = j
            numgb(i+1) = n1
         end do
      end do

c group elements(volume, area, moment of inertia)
      do 100 ie = 1, nei
         nsl = nsf(ie)
         if(icnl.eq.0) then
c     procedure for 3-d problems
            if(nsl.eq.8) then
c              do 99 lei = 1, 8
c                 numgb(lei) = nelm((ie-1)*ns+lei)
c99            continue
c
c renumber local element nodes to conform to convention
c
c gaz don't reverse this is already done (gaz 7-15-09)
            do lei = 1, 8 
              numgb(lei) = nelm((ie-1)*ns+lei)
            enddo
c            do lei = 5, 8 
c              numgb(lei) = nelm((ie-1)*ns+lei-4)
c            enddo
c     3-d 8 node brick
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)
               z12 = cord(numgb(1), 3)-cord(numgb(2), 3)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)
               z14 = cord(numgb(1), 3)-cord(numgb(4), 3)

               x23 = cord(numgb(2), 1)-cord(numgb(3), 1)
               y23 = cord(numgb(2), 2)-cord(numgb(3), 2)
               z23 = cord(numgb(2), 3)-cord(numgb(3), 3)

               x34 = cord(numgb(3), 1)-cord(numgb(4), 1)
               y34 = cord(numgb(3), 2)-cord(numgb(4), 2)
               z34 = cord(numgb(3), 3)-cord(numgb(4), 3)

               x56 = cord(numgb(5), 1)-cord(numgb(6), 1)
               y56 = cord(numgb(5), 2)-cord(numgb(6), 2)
               z56 = cord(numgb(5), 3)-cord(numgb(6), 3)

               x58 = cord(numgb(5), 1)-cord(numgb(8), 1)
               y58 = cord(numgb(5), 2)-cord(numgb(8), 2)
               z58 = cord(numgb(5), 3)-cord(numgb(8), 3)

               x67 = cord(numgb(6), 1)-cord(numgb(7), 1)
               y67 = cord(numgb(6), 2)-cord(numgb(7), 2)
               z67 = cord(numgb(6), 3)-cord(numgb(7), 3)

               x78 = cord(numgb(7), 1)-cord(numgb(8), 1)
               y78 = cord(numgb(7), 2)-cord(numgb(8), 2)
               z78 = cord(numgb(7), 3)-cord(numgb(8), 3)

               x15 = cord(numgb(1), 1)-cord(numgb(5), 1)
               y15 = cord(numgb(1), 2)-cord(numgb(5), 2)
               z15 = cord(numgb(1), 3)-cord(numgb(5), 3)

               x37 = cord(numgb(3), 1)-cord(numgb(7), 1)
               y37 = cord(numgb(3), 2)-cord(numgb(7), 2)
               z37 = cord(numgb(3), 3)-cord(numgb(7), 3)
c     
c     calculate triple products(volumes)
c     
               tp1=tp(-x15,-y15,-z15,-x12,-y12,-z12,-x14,-y14,-z14)
               tp2=tp(-x37,-y37,-z37,-x34,-y34,-z34,x23,y23,z23)
               tp3=tp(-x15,-y15,-z15,-x56,-y56,-z56,-x58,-y58,-z58)
               tp4=tp(-x37,-y37,-z37,-x78,-y78,-z78,x67,y67,z67)
               volest=(tp1+tp2+tp3+tp4)/4.0
               vole(ie)=volest
c     
c     check orthogonality if nonstress solution
c     
               iorth(ie) = 0
            endif
c     
            if(nsl.eq.6) then
               do lei = 1, 6
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     procedure for 6 node prisms
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)
               z12 = cord(numgb(1), 3)-cord(numgb(2), 3)

               x13 = cord(numgb(1), 1)-cord(numgb(3), 1)
               y13 = cord(numgb(1), 2)-cord(numgb(3), 2)
               z13 = cord(numgb(1), 3)-cord(numgb(3), 3)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)
               z14 = cord(numgb(1), 3)-cord(numgb(4), 3)

               x45 = cord(numgb(4), 1)-cord(numgb(5), 1)
               y45 = cord(numgb(4), 2)-cord(numgb(5), 2)
               z45 = cord(numgb(4), 3)-cord(numgb(5), 3)

               x46 = cord(numgb(4), 1)-cord(numgb(6), 1)
               y46 = cord(numgb(4), 2)-cord(numgb(6), 2)
               z46 = cord(numgb(4), 3)-cord(numgb(6), 3)

c     calculate triple products(volumes)
               tp1=tp(-x14,-y14,-z14,-x12,-y12,-z12,-x13,-y13,-z13)
               tp2=tp(-x14,-y14,-z14,-x45,-y45,-z45,-x46,-y46,-z46)
               volest = abs((tp1+tp2)/4.0)
               vole(ie) = volest
c     prisms always are non orthogonal
               iorth(ie) = 0
            endif
            if(nsl.eq.3) then
               do  lei = 1, 3
                  ned = nelm((ie-1)*ns+lei)
                  dm(lei, 1) = cord(ned, 1)
                  dm(lei, 2) = cord(ned, 2)
                  dm(lei, 3) = cord(ned, 3)
               enddo
c     procedure for 3 node triangles in 3-d
c     calculate volume
               call area2d_tri(1,dm(1,1),dm(1,2),dm(1,3),
     &              dm(2,1),dm(2,2),dm(2,3),dm(3,1),dm(3,2),dm(3,3),    
     &              volest,x12,x13,x14,y12,y13,y14)                
               vole(ie) = volest
c     triangles always are non orthogonal
               iorth(ie) = 0
            endif
            if(nsl.eq.4) then
               do lei = 1, 4
                  ned = nelm((ie-1)*ns+lei)
                  dm(lei, 1) = 1.0
                  dm(lei, 2) = cord(ned, 1)
                  dm(lei, 3) = cord(ned, 2)
                  dm(lei, 4) = cord(ned, 3)
               enddo
c     procedure for 4 node tetrahedrals
c     calculate volume
               call determ(det, dm, 4)
               volest = 0.166667*det
               vole(ie) = volest
c     tetrahedrals always are non orthogonal
               iorth(ie) = 0
            endif
         endif
c     
c     procedure for 2-d elements
         if(icnl.ne.0) then
c     4-node quadrilaterals
            if(nsl.eq.4) then
               do lei = 1, 4
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)

               x23 = cord(numgb(2), 1)-cord(numgb(3), 1)
               y23 = cord(numgb(2), 2)-cord(numgb(3), 2)

               x34 = cord(numgb(3), 1)-cord(numgb(4), 1)
               y34 = cord(numgb(3), 2)-cord(numgb(4), 2)

               tp1 = vp(-x12, -y12, -x14, -y14)
               tp2 = vp(-x34, -y34, x23, y23)
               volest = (tp1+tp2)/2.0
               vole(ie) = volest
               iorth(ie) = 0
            endif
            if(nsl.eq.3) then
c     procedure for triangles
               do lei = 1, 3
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)

               x13 = cord(numgb(1), 1)-cord(numgb(3), 1)
               y13 = cord(numgb(1), 2)-cord(numgb(3), 2)

               volest = 0.5*(x12*y13-x13*y12)
               vole(ie) = volest
c     
c     triangles are always nonorthogonal
               iorth(ie) = 0
            endif
         endif
 100  continue
      deallocate(narr,nsf,iorth,nf,iplace)
c
c     ****************** begin anisotropy coding ********************* 
c
c
c     code for anisotropoy calculations
c     this is for 8-node brick elements(icnl=8,ns=8)
c     defined counterclockwise bottom to top
c     since FEHM has non-standard (top to bottom)
c     local numbers are changed for this section
c
c rotate principal axes to x,y,z coordiates
c theta3 is the dip angle
c theta2 is another angle
c theta1 is another angle
c
c here the pnx,pny,pnz represent the principal values of
c permeability in the directions aligned with the stratigraphy
c
      if(ianpe.ne.3) then
      do i = 1,n0
       theta1 = anxy(i)*deg2rad
       theta2 = anxz(i)*deg2rad
       theta3 = anyz(i)*deg2rad
       alen1 = sqrt(1.d00 + cos(theta1 - theta2)**2
     &        *tan(theta3)**2)
       alen2 = sqrt(1.d00 + sin(theta1 - theta2)**2
     &         *tan(theta3)**2)
       a11=cos(theta1)/alen1
       a12=sin(theta1)/alen1
       a13=-cos(theta1-theta2)*tan(theta3)/alen1
       a21=-sin(theta1)/alen2
       a22=cos(theta1)/alen2
       a23=sin(theta1-theta2)*tan(theta3)/alen2
       a31=sin(theta3)*cos(theta2)
       a32=sin(theta3)*sin(theta2)
       a33=cos(theta3)
       ak11 = pnx(i)
       ak22 = pny(i)
       ak33 = pnz(i)
       pnx(i) =ak11*a11*a11+ak22*a21*a21+ak33*a31*a31
       pny(i) =ak11*a12*a12+ak22*a22*a22+ak33*a32*a32
       pnz(i) =ak11*a13*a13+ak22*a23*a23+ak33*a33*a33
       anxy(i)=ak11*a11*a12+ak22*a21*a22+ak33*a31*a32
       anxz(i)=ak11*a11*a13+ak22*a21*a23+ak33*a31*a33
       anyz(i)=ak11*a12*a13+ak22*a22*a23+ak33*a32*a33
c
      enddo
      endif
c
c now permeabilities are transformed
c
      max_elem = 8  
c
c  there are only 4 +x faces per hex
c  the 8 below should therefore be reduced to 4
c
      allocate (sx_temp_x(nei,4,12))
      allocate (sx_temp_y(nei,4,12))
      allocate (sx_temp_z(nei,4,12))
      allocate (ncon_elem(n0,max_elem))
      allocate (ncon_pos(n0,max_elem))
      allocate (iplace(neq))
      allocate (iannip(24,24))
      allocate (ianni(24))
      idani = 3*ns
c    zero out volumes of gridblocks, coefficient storage
      sx_temp_x = 0.0
      sx_temp_y = 0.0
      sx_temp_z = 0.0
      sx1 = 0.0
c ncon_adv not allocated yet
c      ncon_adv = 0
      ncon_elem = 0
      iplace = 0
c
      do ie=1,nei
c
      ianni = 0
      a_ani = 0.0d00
c  
c
c   identify nodes in element
c   note: top and bottom layers are reversed
c   don't need to reverse (gaz 7-15-09)
 
        do lei = 1, 8 
          numgb(lei) = nelm((ie-1)*ns+lei)
        enddo
c        do lei = 5, 8 
c          numgb(lei) = nelm((ie-1)*ns+lei-4)
c        enddo
c   define coeficients for pressure difference
        do j=1,12
         xcell(ithree(j))=(cord(numgb(ione(j)),1)
     &         +cord(numgb(itwo(j)),1))/2.
         ycell(ithree(j))=(cord(numgb(ione(j)),2)
     &         +cord(numgb(itwo(j)),2))/2.
         zcell(ithree(j))=(cord(numgb(ione(j)),3)
     &         +cord(numgb(itwo(j)),3))/2.
        enddo
c test element for isotropy:itest_iso=0
c first set element to isotropic
         itest_iso=0
c
c load permeability data into aklem array
c
         do j=1,ns
          aklem(j,1) = pnx(numgb(j))
          aklem(j,2) = anxy(numgb(j))
          aklem(j,3) = anxz(numgb(j))
          aklem(j,4) = pny(numgb(j))
          aklem(j,5) = anyz(numgb(j))
          aklem(j,6) = pnz(numgb(j))
          if(abs(aklem(j,2))+abs(aklem(j,3))+
     &      abs(aklem(j,5)).gt.a_tol) itest_iso=1
         enddo
c
c !!!!!!! always set itest_iso=1 for now gaz 6-15-04
c
         itest_iso=1
c
c   define area for interface points in control volume 
c
      if(itest_iso.ne.0) then
       call area_interface_hex(1,n0,xcell,ycell,zcell,
     &         ione,itwo,ithree,inine,iten,numgb,cord,area)
      else
       call area_interface_hex_iso(1,n0,xcell,ycell,zcell,
     &         ione,itwo,ithree,ifour,ifive,inine,iten,
     &         numgb,cord,area,8,aklem)
      endif
c
c   form coefficients based on anisotropy 
c
      if(itest_iso.ne.0) then
c
c   load coefficients of transform matrix
c   order: (c1,d1,e1).....(c8,d8,e8)            
c        
        do j=1,12   
         a_ani(j,3*ione(j)-2) = cell_dir(ithree(j))* 
     &         (xcell(ithree(j))-cord(numgb(ione(j)),1))
         a_ani(j,3*ione(j)-1) =  cell_dir(ithree(j))*
     &         (ycell(ithree(j))-cord(numgb(ione(j)),2))
         a_ani(j,3*ione(j)) =  cell_dir(ithree(j))*
     &         (zcell(ithree(j))-cord(numgb(ione(j)),3))
         a_ani(j,3*itwo(j)-2) =  -cell_dir(ithree(j))*
     &        (xcell(ithree(j))-cord(numgb(itwo(j)),1))
         a_ani(j,3*itwo(j)-1) =  -cell_dir(ithree(j))*
     &        (ycell(ithree(j))-cord(numgb(itwo(j)),2))
         a_ani(j,3*itwo(j)) =  -cell_dir(ithree(j))*
     &        (zcell(ithree(j))-cord(numgb(itwo(j)),3))
        enddo
c enforce flux continuity
        do j=1,12    
         jj=j+12
         x2= (cord(numgb(ifour(j)),1)-cord(numgb(ifive(j)),1))**2
         y2= (cord(numgb(ifour(j)),2)-cord(numgb(ifive(j)),2))**2
         z2= (cord(numgb(ifour(j)),3)-cord(numgb(ifive(j)),3))**2
         dis = x2+y2+z2
         cosx=sqrt(x2/dis)
         cosy=sqrt(y2/dis)
         cosz=sqrt(z2/dis)
c
         a_ani(jj,3*ifour(j)-2) =  cosx*aklem(ifour(j),1)
     &     +cosy*aklem(ifour(j),2)+cosz*aklem(ifour(j),3)
         a_ani(jj,3*ifour(j)-1) =  cosx*aklem(ifour(j),2)
     &     +cosy*aklem(ifour(j),4)+cosz*aklem(ifour(j),5)
         a_ani(jj,3*ifour(j))   =  cosx*aklem(ifour(j),3)
     &     +cosy*aklem(ifour(j),5)+cosz*aklem(ifour(j),6)
         a_ani(jj,3*ifive(j)-2) =  -(cosx*aklem(ifive(j),1)
     &     +cosy*aklem(ifive(j),2)+cosz*aklem(ifive(j),3))
         a_ani(jj,3*ifive(j)-1) =  -(cosx*aklem(ifive(j),2)
     &     +cosy*aklem(ifive(j),4)+cosz*aklem(ifive(j),5))
         a_ani(jj,3*ifive(j))   =  -(cosx*aklem(ifive(j),3)
     &     +cosy*aklem(ifive(j),5)+cosz*aklem(ifive(j),6))
c
         a_ani_c(jj,3*ifour(j)) = a_ani(jj,3*ifour(j))
         a_ani_c(jj,3*ifour(j)-1) = a_ani(jj,3*ifour(j)-1)
         a_ani_c(jj,3*ifour(j)-2) = a_ani(jj,3*ifour(j)-2)
         a_ani_c(jj,3*ifive(j)) = a_ani(jj,3*ifive(j))
         a_ani_c(jj,3*ifive(j)-1) = a_ani(jj,3*ifive(j)-1)
         a_ani_c(jj,3*ifive(j)-2) = a_ani(jj,3*ifive(j)-2)
c
        enddo

c
c     ****************** form inverse of a_ani ********************* 
c
c     a_anii(i,j)  will contain its inverse
c
            do j=1,idani
               do k=1,idani
                  a_anii(j,k)=0.0
               enddo
               a_anii(j,j)=1.0
            enddo
            call ludcmp0(a_ani,idani,idani,indx,det)
            if(indx(1).lt.0) then
               indx(1)=-i
               write(ierr, *) 'lu failed in anisotropy calcs: stopping'
               stop    
            endif
            do j=1,idani
               call lubksb0(a_ani,idani,idani,indx,a_anii(1,j))
            enddo
c find non-zeros in a_anii
            do i = 1,idani
             ianni(i)=0
             do k = 1,12
              if(abs(a_anii(i,k)).gt.a_tol) then
               ianni(i)=ianni(i)+1 
               iannip(i,ianni(i))=k
              endif
             enddo
            enddo
c
c     loop over control volume
c
       do kk=1,ns
        vol_con = 0.125*vole(ie)
        sx1(numgb(kk))= sx1(numgb(kk))+vol_con
c      
c       form connectivity matrix
c       
c need to loop over 1 x interface, 1 y interface, 1 z interface
c iface contains the interfaces for eack kk node
c
          do ik=1,3   
           i2=iface(kk,ik)
           j=abs(i2)
           i3 = ik*i2/j
           jj=j+12
           area_j=area(j)
           termx=a_ani_c(jj,3*kk-2)*area_j
           termy=a_ani_c(jj,3*kk-1)*area_j
           termz=a_ani_c(jj,3*kk)*area_j
c     check for other node on face
c     only "+" faces should be necessary
c       changed to local variables 6-2-04
           if(i3.eq.1) then
            kki = kkx(kk)
            if(abs(termx).gt.a_tol) then
             irowc=3*kk-2
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_x(ie,kki,ii)=
     &           sx_temp_x(ie,kki,ii)
     &            +termx*a_anii(irowc,je)
             enddo
            endif
           if(abs(termy).gt.a_tol) then
             irowc=3*kk-1
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_x(ie,kki,ii)=
     &           sx_temp_x(ie,kki,ii)
     &            +termy*a_anii(irowc,je)
             enddo
            endif
           if(abs(termz).gt.a_tol) then
             irowc=3*kk
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_x(ie,kki,ii)=
     &           sx_temp_x(ie,kki,ii)
     &            +termy*a_anii(irowc,je)
             enddo
            endif
           else if(i3.eq.2) then
c       d term (y) (identify row in a_anii)
            kki = kky(kk)
            if(abs(termx).gt.a_tol) then
             irowc=3*kk-2
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_y(ie,kki,ii)=
     &           sx_temp_y(ie,kki,ii)
     &            +termx*a_anii(irowc,je)
             enddo
            endif
           if(abs(termy).gt.a_tol) then
             irowc=3*kk-1
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_y(ie,kki,ii)=
     &           sx_temp_y(ie,kki,ii)
     &            +termy*a_anii(irowc,je)
             enddo
            endif
           if(abs(termz).gt.a_tol) then
             irowc=3*kk
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_y(ie,kki,ii)=
     &           sx_temp_y(ie,kki,ii)
     &            +termz*a_anii(irowc,je)
             enddo
            endif
           else if(i3.eq.3) then
c       e term (z) (identify row in a_anii)
            kki = kkz(kk)
             if(abs(termx).gt.a_tol) then
             irowc=3*kk-2
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_x(ie,kki,ii)=
     &           sx_temp_x(ie,kki,ii)
     &            +termx*a_anii(irowc,je)
             enddo
            endif
           if(abs(termy).gt.a_tol) then
             irowc=3*kk-1
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_z(ie,kki,ii)=
     &           sx_temp_z(ie,kki,ii)
     &            +termy*a_anii(irowc,je)
             enddo
            endif
           if(abs(termz).gt.a_tol) then
             irowc=3*kk
             do k=1,ianni(irowc)
               je=iannip(irowc,k)
                 ii = ithree(je)
                 sx_temp_z(ie,kki,ii)=
     &           sx_temp_z(ie,kki,ii)
     &            +termz*a_anii(irowc,je)
             enddo
            endif

           endif
          enddo
         enddo
c
c summary of face information 
c
c    position 1
c contribution x+(interface 1),y+(interface(5),z+(interface 9)
c    x+: q=area*(a_ani(1,1)*a_anii(1,1)*(p2-p1)
c    position 2
c contribution x-(interface 1),y+(interface(6),z+(interface 10)
c    position 3
c contribution x-(interface 2),y-(interface(6),z+(interface 11)
c    position 4
c contribution x+(interface 2),y-(interface(5),z+(interface 12)
c    position 5
c contribution x+(interface 4),y+(interface(8),z-(interface 9)
c    position 6
c contribution x-(interface 4),y+(interface(7),z-(interface 10)
c    position 7
c contribution x-(interface 3),y-(interface(7),z-(interface 11)
c    position 8
c contribution x+(interface 3),y-(interface(8),z-(interface 12)
        
        endif
c  calculate node-element information
        do  ik = 1, 8
         kb = numgb(ik)
         iplace(kb) = iplace(kb) +1 
         ncon_elem(kb,iplace(kb)) = ie 
         ncon_pos(kb,iplace(kb)) = ik 
        enddo
          
c     ****************** end do loop on elements ******************** 
        enddo
c       close(66)
c     
c     load connectivity array
c
c   gaz 6-7-04 new formulation to save memory
c
      deallocate(iannip,ianni,vole)
      neqp1 = neq+1
      ncon_space = 18*neq
      ncon_space_sx = 18*neq
      allocate(ncon_x1(ncon_space))
      allocate(ncon_x2(ncon_space))
      allocate(ncon_y1(ncon_space))
      allocate(ncon_y2(ncon_space))
      allocate(ncon_z1(ncon_space))
      allocate(ncon_z2(ncon_space))
      allocate(sx_x(ncon_space_sx))
      allocate(sx_y(ncon_space_sx))
      allocate(sx_z(ncon_space_sx))
      allocate(ncon_adv(n0,3))
      ncon_adv = 0

      ncon_x1 = 0
      ncon_y1 = 0
      ncon_z1 = 0
      ncon_x2 = 0
      ncon_y2 = 0
      ncon_z2 = 0

      allocate(icxani(n0+1))
      allocate(icyani(n0+1))
      allocate(iczani(n0+1))

      icx = 0
      icxani(1) = icx
      icy = 00
      icyani(1) = icy
      icz = 0
      iczani(1) = icz


      nedge = 12

      do i = 1,neq
         do jj = 1,max_elem
            ie = ncon_elem(i,jj)
            if(ie.eq.0) go to 700
c identify global and local element numbers
c now done elsewhere gaz 7-14-2009
            do lei = 1, 8 
               numgb(lei) = nelm((ie-1)*ns+lei)
            enddo
c            do lei = 5, 8 
c               numgb(lei) = nelm((ie-1)*ns+lei-4)
c            enddo
            kk = ncon_pos(i,jj)
            do ik = 1,3
c loop on three faces associated with element position
c check for consistency !! 
c  this is coding from last section
c  can possibly remove ncon_adv
               i3 = kb_face_local(kk,ik)      
               if( i3.gt.0) then
                  kb_node = numgb(i3)
                  ncon_adv(i,ik) = kb_node
               else
                  go to 800
               endif
               do je = 1, nedge
c loop on edges of element, identify global nodes
                  ii = ithree(je)
                  kb1 = numgb(ione(je))
                  kb2 = numgb(itwo(je))
                  if(ik.eq.1) then
c x face
                     termx = sx_temp_x(ie,kkx(kk),ii)
                     if(abs(termx).gt.a_tol) then
                        icx = icx +1
                        ncon_x1(icx) = kb1 
                        ncon_x2(icx) = kb2
                        sx_x(icx) = termx
                     endif
                  else if(ik.eq.2) then
c y face
                     termy = sx_temp_y(ie,kky(kk),ii)
                     if(abs(termy).gt.a_tol) then
                        icy = icy +1
                        ncon_y1(icy) = kb1 
                        ncon_y2(icy) = kb2
                        sx_y(icy) = termy
                     endif
                  else  
c z face
                     termz = sx_temp_z(ie,kkz(kk),ii)
                     if(abs(termz).gt.a_tol) then
                        icz = icz +1
                        ncon_z1(icz) = kb1 
                        ncon_z2(icz) = kb2
                        sx_z(icz) = termz
                     endif  
                  endif
               enddo
 800           continue
            enddo
         enddo
 700     continue
         icxani(i+1) = icx
         icyani(i+1) = icy
         iczani(i+1) = icz
      enddo
      deallocate(sx_temp_x,sx_temp_y,sx_temp_z)
      deallocate(iplace,ncon_elem,ncon_pos)
      if(iani_chk) then
c  assemble overall connectivity array from ncon
      if (iptty .ne. 0) write(iptty,*)
      if (iout .ne. 0) write(iout,*)
      sumu = 0
      do i = 1,icx
         sumu = sumu + sx_x(i)
      enddo
      if (iptty .ne. 0) write(iptty,*) 'sum of x coeffs = ', sumu
      if (iout .ne. 0) write(iout,*) 'sum of x coeffs = ', sumu
      sumu =0.0
      do i = 1,icy
         sumu = sumu + sx_y(i)
      enddo
      if (iptty .ne. 0) write(iptty,*) 'sum of y coeffs = ', sumu
      if (iout .ne. 0) write(iout,*) 'sum of y coeffs = ', sumu
      sumu =0.0
      do i = 1,icz
         sumu = sumu + sx_z(i)
      enddo
      if (iptty .ne. 0) then 
         write(iptty,*) 'sum of z coeffs = ', sumu
         write(iptty,*)
      end if
      if (iout .ne. 0) then
         write(iout,*) 'sum of z coeffs = ', sumu
         write(iout,*)
      end if
      endif
c
c  resize arrays
c
c   x face
      allocate(idum(icx))
      idum(1:icx) = ncon_x1(1:icx)
      deallocate(ncon_x1)
      allocate(ncon_x1(icx))
      ncon_x1 = idum
      deallocate(idum)
      allocate(idum(icx))
      idum(1:icx) = ncon_x2(1:icx)
      deallocate(ncon_x2)
      allocate(ncon_x2(icx))
      ncon_x2 = idum
      deallocate(idum)
      allocate(dum(icx))
      dum(1:icx) = sx_x(1:icx)
      deallocate(sx_x)
      allocate(sx_x(icx))
      sx_x = dum
      deallocate(dum)
c   y face
      allocate(idum(icy))
      idum(1:icy) = ncon_y1(1:icy)
      deallocate(ncon_y1)
      allocate(ncon_y1(icy))
      ncon_y1 = idum
      deallocate(idum)
      allocate(idum(icy))
      idum(1:icy) = ncon_y2(1:icy)
      deallocate(ncon_y2)
      allocate(ncon_y2(icy))
      ncon_y2 = idum
      deallocate(idum)
      allocate(dum(icy))
      dum(1:icy) = sx_y(1:icy)
      deallocate(sx_y)
      allocate(sx_y(icy))
      sx_y = dum
      deallocate(dum)
c   z face
      allocate(idum(icz))
      idum(1:icz) = ncon_z1(1:icz)
      deallocate(ncon_z1)
      allocate(ncon_z1(icz))
      ncon_z1 = idum
      deallocate(idum)
      allocate(idum(icz))
      idum(1:icz) = ncon_z2(1:icz)
      deallocate(ncon_z2)
      allocate(ncon_z2(icz))
      ncon_z2 = idum
      deallocate(idum)
      allocate(dum(icz))
      dum(1:icz) = sx_z(1:icz)
      deallocate(sx_z)
      allocate(sx_z(icz))
      sx_z = dum
      deallocate(dum)

c
      allocate (nelm_temp(n0,27))
      nelm_temp = 0
      do ie = 1,nei
       do j1 = 1,8
          kb1 = nelm((ie-1)*ns+j1)
        do j2 = 1,8
          kb2 = nelm((ie-1)*ns+j2)
          ipos1 = ipos(j1,j2)
          ipos2 = ipos(j2,j1)
          nelm_temp(kb1,ipos1) = kb2
          nelm_temp(kb2,ipos2) = kb1
        enddo
       enddo
      enddo
c
c  should set the  max size on nelm here
c     
      deallocate(nelm)
      j = 0
      do i = 1,neq
       do jj = 1,27
        if(nelm_temp(i,jj).ne.0) then
         j = j + 1
        endif
       enddo
      enddo
      ncon_space = j + neqp1
      allocate(nelm(ncon_space))
      allocate(idum(n0))
      ic = neqp1
      nelm(1) = ic
      do i = 1,neq
       minkb = i
       maxkb = i
       idum = 0
       i1 = icxani(i)+1
       i2 = icxani(i+1)
       do jj=i1,i2
        kb1 = ncon_x1(jj)
        kb2 = ncon_x2(jj)
        minkb = min(kb1,kb2,minkb)
        maxkb = max(kb1,kb2,maxkb)
        idum(kb1)= 1
        idum(kb2)= 2
       enddo
       i1 = icyani(i)+1
       i2 = icyani(i+1)
       do jj=i1,i2
        kb1 = ncon_y1(jj)
        kb2 = ncon_y2(jj)
        minkb = min(kb1,kb2,minkb)
        maxkb = max(kb1,kb2,maxkb)
        idum(kb1)= 3
        idum(kb2)= 4
       enddo
       i1 = iczani(i)+1
       i2 = iczani(i+1)
       do jj=i1,i2
        kb1 = ncon_z1(jj)
        kb2 = ncon_z2(jj)
        minkb = min(kb1,kb2,minkb)
        maxkb = max(kb1,kb2,maxkb)
        idum(kb1)= 5
        idum(kb2)= 6
       enddo
c   now the - connections     
       do kj = 1,27
        kb = nelm_temp(i,kj)
        if(kb.ne.0) then
        if(ncon_adv(kb,1).eq.i) then
         i1 = icxani(kb)+1
         i2 = icxani(kb+1)
         do jj=i1,i2
          kb1 = ncon_x1(jj)
          kb2 = ncon_x2(jj)
          minkb = min(kb1,kb2,minkb)
          maxkb = max(kb1,kb2,maxkb)
          idum(kb1)= -1
          idum(kb2)= -2
         enddo         
        endif
        if(ncon_adv(kb,2).eq.i) then
         i1 = icyani(kb)+1
         i2 = icyani(kb+1)
         do jj=i1,i2
          kb1 = ncon_y1(jj)
          kb2 = ncon_y2(jj)
          minkb = min(kb1,kb2,minkb)
          maxkb = max(kb1,kb2,maxkb)
          idum(kb1)= -3
          idum(kb2)= -4
         enddo         
        endif
        if(ncon_adv(kb,3).eq.i) then
         i1 = iczani(kb)+1
         i2 = iczani(kb+1)
         do jj=i1,i2
          kb1 = ncon_z1(jj)
          kb2 = ncon_z2(jj)
          minkb = min(kb1,kb2,minkb)
          maxkb = max(kb1,kb2,maxkb)
          idum(kb1)= -5
          idum(kb2)= -6
         enddo         
        endif
        endif
       enddo
        do jj = minkb,maxkb
         if(idum(jj).ne.0) then
          ic = ic +1
          nelm(ic) = jj
          if(jj.eq.i) nelmdg(i) = ic
         endif
        enddo
       nelm(i+1) = ic
      enddo
      deallocate(nelm_temp,idum)
      allocate(idum(ic))
      idum(1:ic) = nelm(1:ic)
      deallocate(nelm)
      allocate(nelm(ic)) 
      nelm = idum
      deallocate(idum)
c     ****************** end anisotropy coding ********************* 
c     
c     set permeabilities to 1 (incorporated into coefficient array)
      do i=1,neq
       pnx(i) = 1.0
       pny(i) = 1.0
       pnz(i) = 1.0
      enddo
c
c     adjust space for connectivity matrix
c     
c     
c     call md_nodes to break connections for those nodes
c     md_nodes and connection breaks need to be done with 
c     sx_x,sx_y,sx_z,etc.
c     
      call md_nodes(2,0,0)                            
c    
c  we must figure these out 
c      deallocate(ncon_temp_1,ncon_temp_2)
c      deallocate(ncon_face,sx_temp,iannip)
c      deallocate(ianni,idum1,idum2,dum)   
c     
c     insure isotropy  in areas (done above)
c     
c GAZ 02-08-01
c setting neq back 
c
      if(gdpm_flag.ne.0 .or. nriver.ne.0) then
         neq = neq_total
         n0=n0_save
      endif
      if (iout .ne. 0) write(iout, 101)
      if (iptty .ne. 0) write(iptty, 101)
 101  format (' >>>>>> Finished FE coef. Anisotropic calcs <<<<<< ')
      return
      end

      subroutine area_interface_hex(iflg,n0,xcell,ycell,zcell,
     &     ione,itwo,ithree,inine,iten,numgb,cord,area)
c
c calculate the internodal area within 8-node brick
c
      implicit none 
      integer i,iflg,inine(*),iten(*)
      integer j,n0,ione(*),itwo(*),ithree(*),numgb(*)
      real*8 xcell(*),ycell(*),zcell(*),cord(n0,*)
      real*8 dis1,dis2,area(*)
      if(iflg.eq.1) then
c assume orthogonal area
         do i=1,12
            dis1 = sqrt((xcell(inine(i))-xcell(i))**2 + 
     &           (ycell(inine(i))-ycell(i))**2 +
     &           (zcell(inine(i))-zcell(i))**2) 
            dis2 = sqrt((xcell(iten(i))-xcell(i))**2 + 
     &           (ycell(iten(i))-ycell(i))**2 +
     &           (zcell(iten(i))-zcell(i))**2) 
            area(i) = 0.25*dis1*dis2
         enddo
      endif
      return
      end
      subroutine area_interface_hex_iso(iflg,n0,xcell,ycell,zcell,
     &     ione,itwo,ithree,ifour,ifive,
     &     inine,iten,numgb,cord,area,idem,aklem)
c
c calculate the internodal area within 8-node brick
c calculate the internodal area/dis within 8-node brick
c
      implicit none 
      integer i,iflg,idem,inine(*),iten(*)
      integer j,n0,ione(*),itwo(*),ithree(*),numgb(*)
      integer ifour(*),ifive(*)
      real*8 xcell(*),ycell(*),zcell(*),cord(n0,*)
      real*8 aklem(idem,*)                             
      real*8 dis1,dis2,dis3,dis_inv,dis_tol,area(*)
      real*8 x1,x2,y1,y2,z1,z2        
      real*8 permx,permy,permz,cosx,cosy,cosz           
      parameter (dis_tol=1.d-30)
c
c assume orthogonal area
c also calculate internodal distance
      if(iflg.eq.1) then
         do i=1,12
            dis1 = sqrt((xcell(inine(i))-xcell(i))**2 + 
     &           (ycell(inine(i))-ycell(i))**2 +
     &           (zcell(inine(i))-zcell(i))**2) 
            dis2 = sqrt((xcell(iten(i))-xcell(i))**2 + 
     &           (ycell(iten(i))-ycell(i))**2 +
     &           (zcell(iten(i))-zcell(i))**2) 
            x1=cord(numgb(ifour(i)),1)
            x2=cord(numgb(ifive(i)),1)
            y1=cord(numgb(ifour(i)),2)
            y2=cord(numgb(ifive(i)),2)
            z1=cord(numgb(ifour(i)),3)
            z2=cord(numgb(ifive(i)),3)
            dis3 = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2) 
            cosx = sqrt((x2-x1)**2)/dis3
            cosy = sqrt((y2-y1)**2)/dis3
            cosz = sqrt((z2-z1)**2)/dis3
            if(dis3.lt.dis_tol) dis3=dis_tol
            permx =  2.0d00*aklem(ifour(i),1)*aklem(ifive(i),1)
     &           /(aklem(ifour(i),1)+aklem(ifive(i),1)+dis_tol)
            permy =  2.0d00*aklem(ifour(i),4)*aklem(ifive(i),4)
     &           /(aklem(ifour(i),4)+aklem(ifive(i),4)+dis_tol)
            permz =  2.0d00*aklem(ifour(i),6)*aklem(ifive(i),6)
     &           /(aklem(ifour(i),6)+aklem(ifive(i),6)+dis_tol)
            dis_inv= 1.d00/dis3
            area(i) = 0.25*dis1*dis2*dis_inv*
     &           (cosx*permx + cosy*permy + cosz*permz)
         enddo
      endif
      return
      end

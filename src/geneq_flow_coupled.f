      subroutine geneq_flow_coupled()
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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
!
! Author : Sai Rapaka, S Kelkar
!!    Function that numerically computes the derivatives of 
!!    flow and energy residuals due to changes in displacements and 
! pressure from effective stress

      use comai, only: nei, neq, ns
      use combi, only: nelm
      use comei, only: a
      use davidi, only: nmat
      use comfem

      implicit none
      real*8 u_eps
      real*8 Ri(2), Rjup(2), Rjvp(2), Rjwp(2), Rjpp(2)
      integer node_I, node_J, i, j
      integer iau, i_begin, i_end, iel_begin, iel_end

      integer ibug,jbug,kbug,elbug, edgebug
      real*8 epsilon_perm
      epsilon_perm = 1.e-18

      u_eps = 1.0d-6

      do i=1, neq
        node_I = i
        i_begin = nelm(i) + 1
        i_end = nelm(i+1)

        !! Compute the residuals for node_I without changing
        !! the displacements or pressures of any of its neighbors
        call compute_flow_residual(Ri, node_I, 0, 0.0d0, 0.0d0, 0.0d0
     &       ,0.0d0,0)
c.............................................
c fopr debugging s kelkar oct 3 2011
c        i_begin = NodeElems(node_I)+1
c        i_end   = NodeElems(node_I+1)
c        do ibug=iel_begin, iel_end
c           elbug = NodeElems(ibug)
c           do jbug=1,28
c              edgebug = edgeNum1(elbug, jbug)
c              do kbug=1,3
c                 if((permfactor(edgebug,kbug)-1.0).gt.epsilon_perm) then
c                    write(89,*)elbug,edgebug,permfactor(edgebug,kbug)
c                 endif
c              enddo
c           enddo
c        enddo
c.....................................

        do j=i_begin, i_end
          node_J = nelm(j)

          !! For each neighbor, compute the residuals R_i(u_j + du_j),
          !! R_i(v_j + dv_j) and R_i(w_j + dw_j)
          call compute_flow_residual(Rjup, node_I, node_J, u_eps, 
     &         0.0d0, 0.0d0, 0.0d0, 1)
          call compute_flow_residual(Rjvp, node_I, node_J, 0.0d0, 
     &         u_eps, 0.0d0, 0.0d0, 1)
          call compute_flow_residual(Rjwp, node_I, node_J, 0.0d0, 
     &         0.0d0, u_eps, 0.0d0, 1)
          call compute_flow_residual(Rjpp, node_I, node_J, 0.0d0, 
     &         0.0d0, 0.0d0, u_eps, 2)

          !! Numerically compute the derivative dR_i/du_j and so on
          iau = j - (neq + 1)

          a(iau + nmat(1)) = a(iau + nmat(1)) + (Rjpp(1) - Ri(1))/u_eps
          a(iau + nmat(5)) =                    (Rjup(1) - Ri(1))/u_eps
          a(iau + nmat(6)) =                    (Rjvp(1) - Ri(1))/u_eps
          a(iau + nmat(7)) =                    (Rjwp(1) - Ri(1))/u_eps

          a(iau + nmat(3)) = a(iau + nmat(3)) + (Rjpp(2) - Ri(2))/u_eps
          a(iau + nmat(8)) =                    (Rjup(2) - Ri(2))/u_eps
          a(iau + nmat(9)) =                    (Rjvp(2) - Ri(2))/u_eps
          a(iau + nmat(10))=                    (Rjwp(2) - Ri(2))/u_eps
        enddo
      enddo

      end subroutine geneq_flow_coupled

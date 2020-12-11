#ifdef DDFLAG
      	module ddmesh
	
c(SILO)	Must dynamically allocate these to fucntion for complex geometry(more than 1 type of polyhedron)
		use comdti, only : n0
		use comai
		use combi

		PetscInt dd_neq
		Mat A_mat, BP_mat
		PetscInt dd_neq_start, dd_neq_end
		PetscInt ii, jj, kk, icv, ic1_tmp, ic2_tmp, ierr_3
		logical :: x_write = .false., y_write = .false., z_write = .false.
		
		Vec connect
		Vec cordinates
		
		
	      	     						
		
		
   	contains
     	
   		!Function to assign mesh vars to petsc vars
     	subroutine mesh_assign(nelm, iocord) 
     		
     		call VecGetOwnershipRange(neq, 1, dd_neq, ierr_3)
     		
     	   cordinates = iocord
           if (i .ne. neq_primary) backspace il
			do i = 1, neq_read
				read(il,*)
			end do
			icv = 1
			do i = 1, nei_in
				read (il,*) ic1,ic2,char_type,(nelm(j), j=1,ns_in)
					do kk = 1, ns_in						
						connect(icv) = nelm(kk)
						icv = icv + 1					
					end do
				end if
			end do
			close (il)
		end if
     	

     	end subroutine mesh_assign          	   	
		
      end module ddmesh
#endif
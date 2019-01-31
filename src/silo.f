      module silo
	
c(SILO)	Must dynamically allocate these to fucntion for complex geometry(more than 1 type of polyhedron)
	    use comdti, only : n0
	    
			
		character*50 :: silo_meshtype = ''
		character (len=90) :: silofile
		integer silo_dbfile, silo_ierr, silo_err, ndims, nnodes
		integer :: silo_cycl = 1
		integer ii, jj, kk, icv, ic1_tmp, ic2_tmp
		integer scordl, optlistid, optlistid1            
		integer NSHAPETYPES
		parameter (NSHAPETYPES = 1)		                         
		integer shapesize(NSHAPETYPES), shapecounts(NSHAPETYPES), shapetype(NSHAPETYPES)
		real*8 :: silo_time = 0
		logical :: x_write = .false., y_write = .false., z_write = .false.
		
		integer, allocatable :: silo_cv(:)
		

		character(:), allocatable :: silo_conname
		character(:), allocatable :: silo_gdkmname
		character(:), allocatable :: silo_unit
      	integer silo_varname_len
      	integer silo_unit_len
      	      
      	
      	
   	contains
     	
   		!Function to write variables to the silo mesh
     	subroutine silo_put_unstruc2(silo_name, silo_units, silo_array) 
			use comai, only : neq_primary, gdkm_flag, neq_gdkm
			use comdti, only : n0
			use combi, only : igdpm
			use, intrinsic :: IEEE_ARITHMETIC
			use, intrinsic :: iso_fortran_env
			
			implicit none
			
			include "silo.inc"		
			
			integer i
			real*8 :: silo_array(n0)
			real*8 :: silo_gd(neq_primary)
			character(len=*) :: silo_name
			character(len=*) :: silo_units
			
			! Set GDKM array to NaN for visual cleanliness
			silo_gd = IEEE_VALUE(silo_gd,IEEE_QUIET_NAN)
			
			! Set sizes of strings for silo funciton (it needs the length)
			silo_unit_len = len(silo_units)
			silo_varname_len = len(silo_name)
			
			silo_err = dbaddcopt(optlistid, DBOPT_UNITS, silo_units, silo_unit_len)
					
			silo_ierr = dbputuv1(silo_dbfile, silo_name, silo_varname_len, "pointmesh", 9, silo_array, neq_primary,
     &		DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		
     		
     		! Set silo_gd values only for nodes at which gdkm is active
     		if(gdkm_flag .eq. 1) then     			
     			do i = 1, neq_primary
     				if(igdpm(i) .ne. 0) then
     					silo_gd(i) = silo_array(i)
     				end if
     			end do     			
     		end if
     		
     		! Write gdkm values to mesh with gdkm appended to name
     		if(gdkm_flag .eq. 1) then
     			silo_gdkmname = trim(silo_name) // "_gdkm"
     			silo_varname_len = len(silo_gdkmname)
     			print *, " NOTE: gdkm detected for ", silo_gdkmname
     			
				silo_ierr = dbputuv1(silo_dbfile, silo_gdkmname, silo_varname_len, "pointmesh", 9, silo_gd, neq_primary,
     &			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)     			
     		end if     		

     	end subroutine silo_put_unstruc2          	   	
		
      end module silo

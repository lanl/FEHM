      module silo
	
c(SILO)	Must dynamically allocate these to fucntion for complex geometry(more than 1 type of polyhedron)
	    use comdti, only : n0
			
		character*50 :: silo_meshtype = ''
		character (len=90) :: silofile
		integer silo_dbfile, silo_ierr, silo_err, ndims, nnodes
		integer :: silo_cycl = 0
		integer ii, jj, kk, icv, ic1_tmp, ic2_tmp
		integer scordl, optlistid, optlistid1            
		integer NSHAPETYPES
		parameter (NSHAPETYPES = 1)		                         
		integer shapesize(NSHAPETYPES), shapecounts(NSHAPETYPES), shapetype(NSHAPETYPES)
		real :: silo_time = 0		
		logical :: x_write = .false., y_write = .false., z_write = .false.

		real, allocatable :: x_cord(:)
		real, allocatable :: y_cord(:)
		real, allocatable :: z_cord(:)
		integer, allocatable :: silo_cv(:)
		
      end module silo

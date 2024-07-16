c-------- program to make simple grid

      program grid

	real*8  y , x

        open(23,file='grid_out')

	write(23,*) 'coor'
        write(23,*) '    ',12

	do n = 1,12, 2

	  x = (dreal(n)-1.0)*0.1
          write(23,999) n,    0.0, x,  0. 
          write(23,999) n+1,  1.0, x,  0.
	 
	end do

	write(23,*)
	write(23,*) 'elem'
	write(23,*) '    ',4,5

	m = 1

        do n = 1,5
	  write(23,998) n,m, m+1, m+3, m+2
	  m = m + 2
	end do

	        write(23,*)
        write(23,*) 'stop'

        close(23)
998     FORMAT('  ',5I8)
999     FORMAT('  ',I5,3G15.8)
        END



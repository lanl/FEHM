
open(11,file='1dgrid.grid')
open(12,file='flow.00002_sca_node.avs')
read(11,*) 
read(11,*) nnode 

open(21,file='ER-12-3.perm')
open(22,file='ER-12-3.rlp')
open(23,file='ER-12-3.rock')
read(21,*)
read(22,*)
read(23,*)

read(12,*) 
read(12,*) 
read(12,*) 
do i=1,nnode
   read(11,*) n, x, y
   read(12,*) n, h, s
   h = (h-0.1)*1.0e6/997.80831/9.8
      read(21,*) itmp, itmp, itmp, perm
      read(22,*) itmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
      read(23,*) itmp, itmp, itmp, tmp7, tmp8, tmp9
   write(13,'(15e14.5)') x, y, h, s, perm, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
enddo

end

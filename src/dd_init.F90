#ifdef DDFLAG
	subroutine dd_init()
	
            call UGridCreateUGDM(unstructured_grid,ugdm,ndof,option)
            
            call UGridPartition(unstructured_grid,option,Dual_mat,is_new, &
                      num_cells_local_new)
                      
                      
            call UGridDMCreateVector(unstructured_grid,ugdm,global_vec,GLOBAL,option)
            
			call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
			          grid%nG2L,grid%nL2G,grid%nG2A,option)
			          
			          

                      
            call UGridDecompose(discretization%grid%unstructured_grid, &
                              option)                      
                              
            call UGridMapSideSet(grid%unstructured_grid, &
                            region%sideset%face_vertices, &
                            region%sideset%nfaces,region%name, &
                            option,region%cell_ids,region%faces)
                            region%num_cells = size(region%cell_ids)                              

                            
            call DivideVariableArray(data_array, unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)
	end                           
#endif

    1.  Load vertex ids in local arrays.
    2.  MPI_Exscan to get global offset for first cell on processor
    3.  Create and fill MPIAdj (adjacency matrix)
    4.  Call MatMeshToCellGraph() to produce the dual graph
    
    5.  Create and set matrix partitioning based on dual produced in 4.
      a.  MatPartitioningCreate()
      b.  MatPartitioningSetAdjacency()
    6.  Call MatPartitioningApply() to produce an IS (is_new) that indicates the
        processor id of each cell.
    7.  Calculate number of elements on each processor with 
        ISPartitioningCount().
        
    8.  Create a vec (elements_natural) of size num_cells_local_new
    13. Scatter the element vector to the new petsc global
      a.  Create scatter context from is_scatter created above
      
    9.  Create a new IS with the new global number of each index (cell) on
        each processor (is_num)
    10. Create a Blocked IS (is_scatter geh:name not clear)) for a new PetscVec  
    
    20. Take duals and convert to petsc ordering
      a.  Count up number of duals and cells and allocate temporary array
      b.  Set values in array based on elements_natural vec
    22. Load into elements_petsc
    23. Count up the number of ghost cells and place in array
    25. Sort ghost cells and assign back to duals based on local ids
    26. Sort vertex ids and assign to vertex ids natural
    
    27. Create vector and scatter vertices to new processors
    28. Create scatter/gather ISs and perform scatter/gather
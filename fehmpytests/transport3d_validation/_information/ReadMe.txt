This is a 3D solute advective-transport benchmark test in MT3D's user guide, problem 7.7, here named as P77.

We simulated transport of conservative solute (no chemical reactions such as bromide or cloride) injected continuously from a point source in a steady unform flow field. The simulation domain is evenly discretized into 8 layers, 15 rows, and 21 columns. A regular grid spacing of 10 m is used for each row/column/layer. The model layer is simulated as a confined layer. The top and bottom of the model layer are at an elevation of 80 m and 0 m, respectively. The injection point at a rate of 0.5 m3/day and observation points are in the 7th layer

FEHM_3d_finemesh: FEHM model runs with a fine mesh discretization (quad mesh)
FEHM_3d_finemesh/tet_mesh_runs: FEHM model runs with a fine mesh discretization (tetrahedron mesh)

MT3D_P77_3d_finemesh: MT3DMS model runs with a fine mesh
MT3D_P77_3d_finemesh_selectoutputs: MT3DMS model runs with a fine mesh and the selected output times, having exactly same as FEHM's output times
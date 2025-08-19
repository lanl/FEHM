Tests for non-darcy flow 2D

This input directory contains a variety of *.in input files and meshes 
used during development of non darcy capabilities in FEHM.
Tests that are stable and verified have fehmn.files in the control folder.
The control files also have appropiate verified outputs in the compare folder.
fehmpytests will only run files in control.

Update work from George SRC_150625
FEHM V3.6.3.2 DATE QA:NDY        08/07/2025    11:18:01

Update has Same tests but improved results from high ndar values.

Rectangle Test Mesh is 2D XY 20x20 with 1 meter spacing
mdat_2d_simple.grid is regular spaced orthogonal mesh
tri_nonortho.fehmn is same mesh smoothed into unstructured mesh

There is a set of Liquid water two phase and gas tests.
The input files differ by setting ndar to OFF for Darcy flow, 
or beta values 1d9, 1d10, 1d11 each with faster velocities.

Files presWAT are compared within tolerance
Run and View 00002_sca_node.vtk contours to compare visually
Velocities written to .out with nodes 1,421,21 on the Y-axis P(MPa)
Highest pres values at node 1 (0,0) and lowest pres at top node 421 (0,20)

dfn mesh is created with dfnWorks with a test for liquid flow 2D in 3D dimensions

Original Test directories created from George's directories:
FEHM V3.6.3   DATE QA:NDY        06/02/2025    14:15:09

LIQUID
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/FEHMPYTEST_LIQ_tests

Tests verified in George's directory:
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/NON_DARCY_FLOW_2D_TEST_LIQ
/OUTPUT_2D_LIQONLY with darcy, darcy with generated restart, ndar models 1d06, 1d09, 1d10

drwxrwsr-- 3 tamiller sft   4096 Feb 18 10:44 input
drwxrwsr-- 2 tamiller sft   4096 Feb 18 10:40 output_1d06
drwxrwsr-- 2 tamiller sft   4096 Feb 18 10:48 output_1d09
drwxrwsr-- 2 tamiller sft   4096 Feb 18 10:34 output_1d10
drwxrwsr-- 2 tamiller sft   4096 Feb 18 10:23 output_darcy
drwxrwsr-- 2 tamiller sft   4096 Feb 14 17:24 output_darcy_rsto

GAS

example input files:
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/FEHMPYTEST_GAS_tests/output_1d10_gonly

gas_ndar_10.in = non_darcy_2D_gas2_gonly10.in 
compare files from:
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/NON_DARCY_FLOW_2D_TEST_GAS/OUTPUT_2D_GASONLY/non_darcy_2D_gas2_gonly_presWAT.hisND10

Top Folder:
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/FEHMPYTEST_GAS_tests

GAS runs verified by George's files in
/project/eesdev/FEHM/from_george/NONDARCY/NONDARCY_PKG_LANL_280125/NON_DARCY_FLOW_2D_TEST_GAS

/OUTPUT_2D_GASONLY has darcy, models 1d06, 1d10 (gonly)
/OUTPUT_2D has darcy but no non-darcy runs (gas2gn)

drwxrwsr-- 2 tamiller sft 4096 Feb 18 17:35 input
drwxrwsr-- 2 tamiller sft 4096 Feb 19 08:16 output_1d06_gonly
drwxrwsr-- 2 tamiller sft 4096 Feb 18 17:17 output_1d10_gas2gn
drwxrwsr-- 2 tamiller sft 4096 Feb 18 17:33 output_1d10_gonly
drwxrwsr-- 2 tamiller sft 4096 Feb 18 16:52 output_darcy_gas2gn
drwxrwsr-- 2 tamiller sft 4096 Feb 18 17:24 output_darcy_gonly
-rw-rw-r-- 1 tamiller sft 1889 Feb 18 17:22 README.txt


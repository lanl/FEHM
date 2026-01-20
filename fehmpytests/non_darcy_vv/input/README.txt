
ND Tests based on /project/eesdev/FEHM/VV_TESTS/V+V_TESTS_ND_050126

Use gaz Windows output files for comparisons.
V+V_TESTS_NDAWH/OUTPUT_NDAWH_gaz/test_ND_AWH_LV_ND12_satr.his
V+V_TESTS_NDAWH/OUTPUT_NDAWH_gaz/test_ND_AWH_LV_NDOFF_satr.his
V+V_TESTS_NDWH/OUTPUT_NDWH_gaz/test_ND_WH_LV_ND12_satr.his
V+V_TESTS_NDW2P/OUTPUT_NDW2P_gaz/test_NDW2Pa_ND12.00002_sca_node.avs

V+V_TESTS_NDWH/OUTPUT_NDWH_gaz/test_ND_WH_LV_NDOFF_satr.his (not used since OFF was not set)
V+V_TESTS_NDW2P/OUTPUT_NDW2P_gaz/test_NDW2Pa_ND00.00002_sca_node.avs (not used, setting unknown)

From: george zyvoloski <gazyvoloski@gmail.com>
Sent: Saturday, January 3, 2026 2:19 PM

I sent test files for ND with AWH, WH, and AW2P. 
They are 1D and designed to test all phases.
They are rather large because I included a folder with EOS tables. 
The attached word document describes the tests.

These should be designated as "truth". 
Each contains the result of two simulations: ND=OFF and ND=1.D12. 
Simply toggle the "off" in the input file:
ndar  off
1 0 0 1.d12

end ndar

For future comparisons, I would use the history files like:
test_ND_AWH_LV_ND12_satr.his for transient comparison
test_ND_AWH_LV_ND12.00002_sca_node.dat for spatial comparison.


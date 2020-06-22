## FEHM-SILO HOW TO BUILD AND USE

This software and all dependencies are open source software available under the BSD-3 or BSD license.
 
### To obtain the latest (passing) version of FEHM-SILO:

git clone https://github.com/lanl/FEHM.git --branch SILO  
cd FEHM/src/  
make  

### To use SILO:
open the fehm.dat file  
change the output method ( 'avs', 'tec', etc ) to 'silo'  
insure 'xyz' , 'geom' , and 'mat , are all in the cont macro  

view the 'material.silo' , 'concen####.silo' , or 'scalar####.silo' files using paraview or visit  

ignore any other .silo files  

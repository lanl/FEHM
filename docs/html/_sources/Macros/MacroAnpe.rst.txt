========
``anpe``
========

Anisotropic permeability input. Adds cross terms to the perm macro. 

The ANPE keyword implements a flux-continuous anisotropic permeability tensor with cross terms. The cross terms can either be input directly or grid rotation angles inputted and the cross terms calculated by FEHM. FEHM implements the method presented by Lee et al 2002. 

The ANPE is incompatible with keywords GDKM, GDPM, DUAL, and DPDP.

* `Group 1 - JA, JB, JC, ANXY, ANXZ, ANYZ 

   - (JA, JB, JC - `are defined here <Macro20058.html>`_)

+----------------+--------+---------+----------------------------------------------------+
| Input Variable | Format | Default | Description                                        |
+================+========+=========+====================================================+
| ANXY           | real   | 1.e-30  | Anisotropic permeability in the xy-direction (m2). |
+----------------+--------+---------+----------------------------------------------------+
| ANXZ           | real   | 1.e-30  | Anisotropic permeability in the xz-direction (m2). |
+----------------+--------+---------+----------------------------------------------------+
| ANYZ           | real   | 1.e-30  | Anisotropic permeability in the yz-direction (m2). |
+----------------+--------+---------+----------------------------------------------------+

#. Lee et al., 2002, Implementation of a Flux-Continuous Finite-Difference Method for Stratigraphic Hexahedron Grids, SPE Journal, Volume 7, Number 3, DOI 10.2118/80117-PA.
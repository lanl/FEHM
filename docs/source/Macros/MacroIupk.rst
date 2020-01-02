========
``iupk``
========

No input is associated with this control statement. This macro enables upwinding, the technique of evaluating the non-linear equation coefficients using the direction of flow relative to the grid block. For example, if flow is moving from grid block j to i, the coefficients for block i, are evaluated at the "upwind" block j. When upwinding is enabled the full transmissibility term will be upwinded (including the intrinsic permeability). Otherwise the fluid and relative permeability part of the transmissibility will be upwinded and the intrinsic permeability will be harmonically averaged.


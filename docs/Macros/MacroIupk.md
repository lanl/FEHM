---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">iupk</span></code><a class="headerlink" href="#iupk" title="Permalink to this headline">¶</a></h1>
<p>No input is associated with this control statement. This macro enables upwinding, the technique of evaluating the non-linear equation coefficients using the direction of flow relative to the grid block. For example, if flow is moving from grid block j to i, the coefficients for block i, are evaluated at the “upwind” block j. When upwinding is enabled the full transmissibility term will be upwinded (including the intrinsic permeability). Otherwise the fluid and relative permeability part of the transmissibility will be upwinded and the intrinsic permeability will be harmonically averaged.</p>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
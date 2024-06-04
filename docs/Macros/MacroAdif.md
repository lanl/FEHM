---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">adif</span></code><a class="headerlink" href="#adif" title="Permalink to this headline">¶</a></h1>
<p>Air-water vapor diffusion. The air-water diffusion equation is given as Equation (21)
of the “Models and Methods Summary” of the FEHM Application (Zyvoloski et al. 1999).</p>
<ul class="simple">
<li>Group 1 - TORT</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="8%" />
<col width="76%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>TORT</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Tortuosity for air-water vapor diffusion.</div>
<div class="line">If TORT &gt; 0, <span class="math notranslate nohighlight">\(\tau\)</span> of eqn 21, otherwise</div>
<div class="line">if TORT &lt; 0, <span class="math notranslate nohighlight">\(abs(\tau \phi S_v)\)</span> of the same equation.</div>
<div class="line">If TORT &gt; 1, water-vapor diffusion coefficient is set equal to the value</div>
<div class="line">specified for the first vapor species defined in the trac macro.</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">adif</span></code>. In this example the tortuosity (<span class="math notranslate nohighlight">\(\tau\)</span>)
for vapor diffusion is specified to be 0.8.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>adif</td>
</tr>
<tr class="row-even"><td>0.8</td>
</tr>
</tbody>
</table>

  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>

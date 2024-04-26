---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">hflx</span></code><a class="headerlink" href="#hflx" title="Permalink to this headline">¶</a></h1>
<p>Heat flow input.</p>
<p>A negative heat flow indicates heat flow into the reservoir.</p>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="7%" />
<col width="8%" />
<col width="70%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>QFLUX</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">If QFLXM = 0, then QFLUX is the heat flow (MW). If QFLXM ≠ 0,
then QFLUX is a temperature (oC) and the heat flow is
calculated according to the formula:</div>
<div class="line"><span class="math notranslate nohighlight">\(Q_H = QFLXM \cdot (T-QFLUX)\)</span> (MW).</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>QFLXM</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">If QFLXM &gt; 0, multiplier for heat flow equation given in QFLUX
description (MW/oC). This must be large for large temperature gradients,
or when a constant temperature must be maintained.</div>
<div class="line">If QFLXM &lt; 0, then QFLUX is interpreted as a fixed saturation and</div>
<div class="line"><span class="math notranslate nohighlight">\(Q_H = ABS(QFLXM) \cdot (S_l - QFLUX)\)</span> (MW).</div>
</div>
</td>
</tr>
</tbody>
</table>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">If the parameter <code class="docutils literal notranslate"><span class="pre">QFLUX</span></code> is set to -999, there will be fixed temperature from restart.</p>
</div>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">hflx</span></code>. In this example, at each node from 401 to 410,
a heat flow of 0.001 MW is being injected into the model.</p>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="19%" />
<col width="11%" />
<col width="30%" />
<col width="19%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>hflx</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>401</td>
<td>410</td>
<td>1</td>
<td>-0.001</td>
<td>0.0</td>
</tr>
</tbody>
</table>
<p>In this example, <code class="docutils literal notranslate"><span class="pre">QFLUX</span></code> is set to -999 (see note above).</p>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="19%" />
<col width="11%" />
<col width="30%" />
<col width="19%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>hflx</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>-1</td>
<td>0</td>
<td>0</td>
<td>-999</td>
<td>1e6</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
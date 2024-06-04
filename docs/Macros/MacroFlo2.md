---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">flo2</span></code><a class="headerlink" href="#flo2" title="Permalink to this headline">¶</a></h1>
<ul class="simple">
<li>Group 1 - JA, JB, JC, JD, SKD, EFLOW, AIPED</li>
</ul>
<p>Multiple lines of input may be used, terminated by a blank line.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="8%" />
<col width="77%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>JA</td>
<td>integer</td>
<td rowspan="4">Indices used to define planes in a 3-D simulation with a regular numbering pattern.</td>
</tr>
<tr class="row-odd"><td>JB</td>
<td>integer</td>
</tr>
<tr class="row-even"><td>JC</td>
<td>integer</td>
</tr>
<tr class="row-odd"><td>JD</td>
<td>integer</td>
</tr>
</tbody>
</table>
<p>The flow rates are defined within the inner loop of the do loops:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span>DO JK = JA, JB
  KL = JK - JA
  DO IJ = JA + KL, JC + KL, JD
      ⋅ ⋅ ⋅
  ENDDO
ENDDO
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">elem</span></code><a class="headerlink" href="#elem" title="Permalink to this headline">¶</a></h1>
<p>Element connectivity data. These data are usually created by a mesh generation program, then cut and copied into the input file or a separate geometry data input file.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last"><code class="docutils literal notranslate"><span class="pre">elem</span></code> is required if macro <code class="docutils literal notranslate"><span class="pre">fdm</span></code> is not used.</p>
</div>
<ul class="simple">
<li>Group 1 - <code class="docutils literal notranslate"><span class="pre">NS</span></code>, <code class="docutils literal notranslate"><span class="pre">NEI</span></code></li>
<li>Group 2 - <code class="docutils literal notranslate"><span class="pre">MB</span></code>, <code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(1)</span></code>, <code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(2)</span></code>, …, <code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(NS)</span></code></li>
</ul>
<p>If <span class="math notranslate nohighlight">\(NS &lt; 0\)</span> then <span class="math notranslate nohighlight">\(ABS(NS)\)</span> is interpreted as the number of nodes per element.
<span class="math notranslate nohighlight">\(NS &lt; 0\)</span> signals the code to make rectangles (or bricks in three dimensions)
a sum of triangles (or tetrahedrals). This provides more stability in nonlinear problems
with a distorted mesh. <a class="reference external" href="Macro39151.html">See Elements available with FEHM in 2-D and 3-D problems showing
nodal numbering convention.</a> shows available element types
and the nodal numbering convention. To end the control section a blank line is entered.</p>
<table border="1" class="docutils">
<colgroup>
<col width="7%" />
<col width="4%" />
<col width="88%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>Input Variable</td>
<td>Format</td>
<td>Description</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">NS</span></code></td>
<td>integer</td>
<td>Number of nodes per element.</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">NEI</span></code></td>
<td>integer</td>
<td>Number of elements</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">MB</span></code></td>
<td>integer</td>
<td>Element number. <code class="docutils literal notranslate"><span class="pre">If</span> <span class="pre">MB</span> <span class="pre">&lt;</span> <span class="pre">0</span></code> then the difference between the absolute value of MB and the previous absolute value of MB is used to generate intermediate values by interpolation in the code.</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(1)</span></code></td>
<td>integer</td>
<td>First node of element MB</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(2)</span></code></td>
<td>integer</td>
<td>Second node of element MB</td>
</tr>
<tr class="row-odd"><td colspan="3">…</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">NELM</span> <span class="pre">(NS)</span></code></td>
<td>integer</td>
<td>Last node of element MB</td>
</tr>
</tbody>
</table>
<p>The following is an example of <strong>elem</strong>. In this example there are 4 nodes per element,
i.e., the elements are 2-dimensional quadrilaterals. There are a total of
117 elements in the model, element number 1 is defined by nodes 15, 16, 2, and 1,
element number 2 is defined by nodes 16, 17, 3 and2, and so on.</p>
<table border="1" class="docutils">
<colgroup>
<col width="23%" />
<col width="19%" />
<col width="19%" />
<col width="19%" />
<col width="19%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>elem</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>4</td>
<td>117</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>15</td>
<td>16</td>
<td>2</td>
<td>1</td>
</tr>
<tr class="row-even"><td>2</td>
<td>16</td>
<td>17</td>
<td>3</td>
<td>2</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>10</td>
<td>24</td>
<td>25</td>
<td>11</td>
<td>10</td>
</tr>
<tr class="row-odd"><td>11</td>
<td>25</td>
<td>26</td>
<td>12</td>
<td>11</td>
</tr>
<tr class="row-even"><td>12</td>
<td>26</td>
<td>27</td>
<td>13</td>
<td>12</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>116</td>
<td>138</td>
<td>139</td>
<td>125</td>
<td>124</td>
</tr>
<tr class="row-odd"><td>117</td>
<td>139</td>
<td>140</td>
<td>126</td>
<td>125</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
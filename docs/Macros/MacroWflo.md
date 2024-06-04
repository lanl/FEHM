---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">wflo</span></code><a class="headerlink" href="#wflo" title="Permalink to this headline">¶</a></h1>
<p>Create a new flow macro to represent boundary conditions on an extracted submodel. Alternate submodel boundary output.</p>
<ul class="simple">
<li>Group 1 - KEYMS1, KEYMS2, KEYMS3, KEYMS3</li>
</ul>
<p>or, if KEYMS4 = ‘type’:</p>
<ul class="simple">
<li>Group 1 - KEYMS1, KEYMS2, KEYMS3, KEYMS4, ITYPSD</li>
</ul>
<p>Enter one line per model defined, terminated with a blank line.</p>
<ul class="simple">
<li>Group 2 - JA, JB, JC, ISUBMD</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="10%" />
<col width="76%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYMS1</td>
<td>character</td>
<td><div class="first last line-block">
<div class="line">Macro SKD type. The letters given in ( ) are sufficient to identify the keyword:</div>
<div class="line">(p)ressure</div>
<div class="line">(h)ead</div>
<div class="line">(f)lux</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KEYMS2</td>
<td>character</td>
<td><div class="first last line-block">
<div class="line">Macro ESKD type. The letters given in ( ) are sufficient to identify the keyword:</div>
<div class="line">(s)aturation</div>
<div class="line">(t)emperature</div>
<div class="line">(e)nthalpy</div>
<div class="line">(w)ater (water only source, output saturation as 1.0)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>KEYMS3</td>
<td>character</td>
<td><div class="first last line-block">
<div class="line">Macro AIPED type:</div>
<div class="line">imph  - Impedance parameter = 1.0</div>
<div class="line">impl - Impedance parameter = 1.0e-4</div>
<div class="line">impn  - Impedance parameter = -1.0</div>
<div class="line">If KEYSM1 = ‘flux’ the impedance parameter = 0.0, otherwise for any other input,
the impedance parameter = 1.0e2.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KEYMS4</td>
<td>character</td>
<td><div class="first last line-block">
<div class="line">Flag to indicate submodel type will be input:</div>
<div class="line">type</div>
</div>
</td>
</tr>
<tr class="row-even"><td>ITYPSD</td>
<td>INTEGER</td>
<td><div class="first last line-block">
<div class="line">Submodel type:</div>
<div class="line"><span class="math notranslate nohighlight">\(ITYPSD = 0\)</span>, generate ‘flow’ macro.</div>
<div class="line"><span class="math notranslate nohighlight">\(ITYPSD \ne 0\)</span>, generate ‘flo3’ macro.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>ISUBMD</td>
<td>integer</td>
<td>Submodel assignment.</td>
</tr>
</tbody>
</table>
<p>The following is an example of the wflo macro. A flow macro will be written for two output zones.
Because, they are specified using two models, the macros will be written to two separate files.
The macro file names will be generated using the root file name (as input or determined from the output
file name) appended with the model number and suffix <code class="docutils literal notranslate"><span class="pre">.wflow</span></code> (e.g. <code class="docutils literal notranslate"><span class="pre">file.0001.wflow</span></code> and <code class="docutils literal notranslate"><span class="pre">file.0002.wflow</span></code>).</p>
<table border="1" class="docutils">
<colgroup>
<col width="29%" />
<col width="24%" />
<col width="29%" />
<col width="19%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>wflo</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>flux</td>
<td>wat</td>
<td>impl</td>
<td>na</td>
</tr>
<tr class="row-odd"><td>flux</td>
<td>wat</td>
<td>impl</td>
<td>na</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>-308</td>
<td>0</td>
<td>0</td>
<td>1</td>
</tr>
<tr class="row-even"><td>-309</td>
<td>0</td>
<td>0</td>
<td>2</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
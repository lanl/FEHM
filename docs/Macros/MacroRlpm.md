---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rlpm</span></code><a class="headerlink" href="#rlpm" title="Permalink to this headline">¶</a></h1>
<p>Relative permeability and capillary pressure model. Several models are available.</p>
<p>Group 1 - KEYWORD “group”, GROUP_NUMBER</p>
<ul class="simple">
<li>The Group 1 KEYWORD “group”, which starts each model sequence.</li>
</ul>
<p>Group 2 - KEYWORD “rlp” PHASE, MODEL_TYPE, RLP_PARAM(i), i = 1, NUMP</p>
<ul class="simple">
<li>At minimum, one entry is made for relative permeability of the wetting phase in the model.  If an entry for a non-wetting phase is not present, Rnw will be calculated using the same model and parameters as Rw.</li>
</ul>
<p>Group 3 (optional) - KEYWORD “cap”, COUPLE, MODEL_TYPE, CAP_PARAM(i), i = 1, NUMP</p>
<ul class="simple">
<li>One entry is made for each phase-couple in the model.  If absent, Cp=0</li>
</ul>
<p>NUMP is the number of parameters needed for the selected model type (see parameter table).</p>
<p>Alternatively, to enter relative permeability and capillary pressure values in a table.</p>
<p>Group 2 -       KEYWORD “table”,  COUPLE</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">For backwards compatability with earlier versions of FEHM Group 2 can be: KEYWORD “table” Table#  PHASE1 PHASE2 COUPLE. Table#,PHASE1,&amp; PHASE2 will be ignored.</p>
</div>
<p>Group 3 -       SATURATION, PHASE1 RELPERM, PHASE2 RELPERM, CAPILLARY PRESSURE</p>
<p>Group 3 is entered multiple times to cover the full range of wetting phase saturation (0. - 1.0).</p>
<p>Table input is terminated with KEYWORD ‘end’ or a blank line. Note input saturation is the saturation of the wetting phase, and PHASE1 is the wetting phase.</p>
<p>-or-</p>
<p>Group 3 -       KEYWORD “file”</p>
<blockquote>
<div>TABLE_FILE</div></blockquote>
<p>Table data will be read from the specified file.</p>
<p>Model input is terminated with KEYWORD ‘end’ or a blank line.</p>
<p>Group 4 -       JA, JB, JC, GROUP_NUMBER</p>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="9%" />
<col width="78%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>group rlp cap</td>
</tr>
<tr class="row-odd"><td>GROUP_NUMBER</td>
<td>integer</td>
<td>Identifier for the current model.</td>
</tr>
<tr class="row-even"><td>PHASE</td>
<td>character</td>
<td>Fluid state: water air co2_liquid co2_sc co2_gas vapor methane_hydrate oil gas</td>
</tr>
<tr class="row-odd"><td>MODEL_TYPE</td>
<td>character</td>
<td>See Tables 1-3</td>
</tr>
<tr class="row-even"><td>RLP_PARAM</td>
<td>real</td>
<td>Input parameters for the specified relative permeability model. See Tables 1-2</td>
</tr>
<tr class="row-odd"><td>NUMP</td>
<td>integer</td>
<td>Number of parameters specified for the model (see Tables 1 and 3).</td>
</tr>
<tr class="row-even"><td>COUPLING</td>
<td>character</td>
<td>Phase couple air/water water/co2_liquid water/co2_gas co2_liquid/co2_gas water/vapor air/vapor</td>
</tr>
<tr class="row-odd"><td>CAP_PARAM</td>
<td>real</td>
<td>Input parameters for the specified capillary pressure model. See Table 3</td>
</tr>
<tr class="row-even"><td>TBLNUM</td>
<td>integer</td>
<td>Table identifier.</td>
</tr>
<tr class="row-odd"><td>PHASE1</td>
<td>character</td>
<td>Fluid state of the wetting phase (see PHASE for phase ID)</td>
</tr>
<tr class="row-even"><td>PHASE2</td>
<td>character</td>
<td>Fluid state of the non-wetting phase (see PHASE for phase ID)</td>
</tr>
</tbody>
</table>
<p>Optional parameters are shown in <strong>bold</strong>.
Table 1. Relative permeability options</p>
<table border="1" class="docutils">
<colgroup>
<col width="57%" />
<col width="20%" />
<col width="7%" />
<col width="4%" />
<col width="12%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Model type</th>
<th class="head">Input Parameter</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>&#160;</td>
<td>1</td>
<td>2</td>
<td>3</td>
<td>4</td>
</tr>
<tr class="row-odd"><td>Linear</td>
<td>Sr</td>
<td>Smax</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Exponential</td>
<td>Sr</td>
<td>Smax</td>
<td>λ</td>
<td><strong>C</strong></td>
</tr>
<tr class="row-odd"><td>Corey</td>
<td>Sr</td>
<td>Smax</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Brooks-Corey</td>
<td>Sr</td>
<td>Smax</td>
<td>λ</td>
<td><strong>D</strong></td>
</tr>
<tr class="row-odd"><td>Van Genuchten (5 model type options, Table 2)</td>
<td>Sr</td>
<td>Smax</td>
<td>m</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>Table 2.  Van Genuchten Relative permeability options</p>
<table border="1" class="docutils">
<colgroup>
<col width="33%" />
<col width="33%" />
<col width="34%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Keyword</th>
<th class="head">Independent variable</th>
<th class="head"><span class="math notranslate nohighlight">\(R_{nw}\)</span></th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><span class="math notranslate nohighlight">\(Vg_1\)</span></td>
<td>S</td>
<td><span class="math notranslate nohighlight">\(1.-R_{w}\)</span></td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(Vg\)</span></td>
<td>S</td>
<td>VG formula</td>
</tr>
<tr class="row-even"><td><span class="math notranslate nohighlight">\(Vg_{1_{cap}}\)</span></td>
<td>Cp</td>
<td><span class="math notranslate nohighlight">\(1.-R_w\)</span></td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(Vg_{cap}\)</span></td>
<td>Cp</td>
<td>Roseangela</td>
</tr>
<tr class="row-even"><td><span class="math notranslate nohighlight">\(Vg_{corey}\)</span></td>
<td>S</td>
<td>Corey formula</td>
</tr>
</tbody>
</table>
<p>For dual porosity problems with Van Genuchten models, an optional line can follow the ‘rlp’ line with the following format:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">KEYWORD</span> <span class="s2">&quot;fracture&quot;</span> <span class="n">Sr</span> <span class="n">Smax</span> <span class="n">m</span> <span class="n">Fracture_k</span> <span class="n">Matrix_k</span> <span class="n">Vf</span> <span class="n">FMIT</span>
</pre></div>
</div>
<p>Note: for backwards compatability with earlier versions of FEHM:</p>
<ol class="arabic simple">
<li>the following input style is accepted for all the Van Genuchten model. <span class="math notranslate nohighlight">\(α_G\)</span> will be ignored.</li>
</ol>
<table border="1" class="docutils">
<colgroup>
<col width="34%" />
<col width="9%" />
<col width="14%" />
<col width="30%" />
<col width="14%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>Van Genuchten</td>
<td>Sr</td>
<td>Smax</td>
<td><span class="math notranslate nohighlight">\(α_G\)</span></td>
<td>n</td>
</tr>
</tbody>
</table>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">KEYWORD</span> <span class="s2">&quot;fracture&quot;</span> <span class="n">Sr</span> <span class="n">Smax</span> <span class="n">m</span> <span class="n">α_G</span> <span class="n">Fracture_k</span> <span class="n">Matrix_k</span> <span class="n">Vf</span> <span class="n">FMIT</span>
</pre></div>
</div>
<ol class="arabic" start="2">
<li><p class="first">the keyword ‘rlp’ can be omitted.</p>
<p>Table 3. Capillary Pressure options</p>
</li>
</ol>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="7%" />
<col width="22%" />
<col width="20%" />
<col width="21%" />
<col width="5%" />
<col width="5%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Model type</th>
<th class="head">&#160;</th>
<th class="head">Input Parameter</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>&#160;</td>
<td>1</td>
<td>2</td>
<td>3</td>
<td>4</td>
<td>5</td>
<td>6</td>
</tr>
<tr class="row-odd"><td>Linear</td>
<td>Sr</td>
<td>Smax</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Linear_for</td>
<td>Cpo</td>
<td>Sco</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>Exponential</td>
<td>Sr</td>
<td>Smax</td>
<td>C</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Brooks-Corey</td>
<td>Sr</td>
<td>Smax</td>
<td>λ</td>
<td><span class="math notranslate nohighlight">\(C_{pe}\)</span></td>
<td>S1</td>
<td>S2</td>
</tr>
<tr class="row-odd"><td>Van Genuchten</td>
<td>Sr</td>
<td>Smax</td>
<td><span class="math notranslate nohighlight">\(α_G\)</span></td>
<td>n</td>
<td>S1</td>
<td>S2</td>
</tr>
</tbody>
</table>
<p>For dual porosity problems with Van Genuchten models, an optional line can follow the ‘cap’ line with the following format:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">KEYWORD</span> <span class="s2">&quot;fracture&quot;</span> <span class="n">Sr</span> <span class="n">Smax</span> <span class="n">m</span>
</pre></div>
</div>
<p>Table 4</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="7%" />
<col width="75%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Parameter</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>Sr</td>
<td>real</td>
<td>Residual (or minimum) wetting-phase saturation.</td>
</tr>
<tr class="row-odd"><td>Smax</td>
<td>real</td>
<td>Maximum wetting-phase saturation.</td>
</tr>
<tr class="row-even"><td>λ</td>
<td>real</td>
<td>Exponent in exponential &amp; Brooks-Corey models</td>
</tr>
<tr class="row-odd"><td>α</td>
<td>real</td>
<td>Inverse of air entry head, <span class="math notranslate nohighlight">\(α_G\)</span> (1/m)</td>
</tr>
<tr class="row-even"><td>n</td>
<td>real</td>
<td>Parameter n from van Genuchten (1980).</td>
</tr>
<tr class="row-odd"><td>C</td>
<td>real</td>
<td>Optional constant in exponential model, =1 if omitted</td>
</tr>
<tr class="row-even"><td>k</td>
<td>real</td>
<td>Optional constant in Brooks-Corey model, (Miller et al., 1998, pg.88) =2 if omitted</td>
</tr>
<tr class="row-odd"><td>Cpo</td>
<td>real</td>
<td>Capillary pressure at zero saturation (MPa)</td>
</tr>
<tr class="row-even"><td>So</td>
<td>real</td>
<td>Saturation when Capillary pressure is zero (-)</td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(C_{pe}\)</span></td>
<td>real</td>
<td>Capillary entry pressure (MPa)</td>
</tr>
<tr class="row-even"><td>Fracture_k</td>
<td>real</td>
<td>Fracture permeability (m2)</td>
</tr>
<tr class="row-odd"><td>Matrix_k</td>
<td>real</td>
<td>Matrix permeability (m2)</td>
</tr>
<tr class="row-even"><td>Vf</td>
<td>real</td>
<td>Fracture volume fraction</td>
</tr>
<tr class="row-odd"><td>FMIT</td>
<td>real</td>
<td>Fracture-matrix interaction term.</td>
</tr>
<tr class="row-even"><td>S1</td>
<td>real</td>
<td>VG fitting parameter</td>
</tr>
<tr class="row-odd"><td>S2</td>
<td>real</td>
<td>VG cutoff saturation</td>
</tr>
</tbody>
</table>
<p>Both Brooks-Corey and Van Genuchten models can be unstable at low saturations.  Parameters S1 and S2 can be used to provide approximations to the curves at low saturations, avoiding numerical instability.</p>
<p>The approximation will replace the B-C or V-G model for saturations less than S2.</p>
<p>If S &lt; S2, then:</p>
<ul class="simple">
<li>If S1=0  (Van Genuchten model only), the curve will be forced to have zero slope at S=0</li>
<li>If S1&gt;0, Cp will vary linearly from Cp(S=S2) to Cp(S=0)=S1*Cp(S=S2) where Cp(S=S2) is capillary pressure at S=S2 calculated according to the B-C or V-G model</li>
<li>If S1&lt;0, Cp will vary linearly from Cp(S=S2) to Cp(S=0) at the slope of the Cp curve at S2 calculated according to the B-C or V-G model</li>
</ul>
<div class="section" id="examples">
<h2>Examples<a class="headerlink" href="#examples" title="Permalink to this headline">¶</a></h2>
<p><strong>Example 1</strong>. Water Rw according to corey, non-wetting phase will be 1-Rw. Cp=0</p>
<table border="1" class="docutils">
<colgroup>
<col width="23%" />
<col width="23%" />
<col width="23%" />
<col width="16%" />
<col width="16%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">rlpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>group</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>rlp</td>
<td>water</td>
<td>corey</td>
<td>0.3</td>
<td>0.1</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>140</td>
<td>1</td>
<td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p><strong>Example 2</strong>. Air/Water problem.  Both water and air will have linear rlp models.  Capillary pressure of air/water will be
according to the linear forsyth model.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="24%" />
<col width="27%" />
<col width="13%" />
<col width="13%" />
<col width="7%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">rlpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>group</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>rlp</td>
<td>water</td>
<td>linear</td>
<td>0.3</td>
<td>1.0</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>rlp</td>
<td>air</td>
<td>linear</td>
<td>0.1</td>
<td>.7</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>cap</td>
<td>air/water</td>
<td>linear_for</td>
<td>93.6</td>
<td><ol class="first last arabic simple" start="100">
<li></li>
</ol>
</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p><strong>Example 3</strong>.  Air/Water problem.  Water (and air, by default) will both have vg rel perm model (see Table 2 for details).
Air/water capillary pressure will be vg.</p>
<table border="1" class="docutils">
<colgroup>
<col width="12%" />
<col width="18%" />
<col width="10%" />
<col width="13%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="7%" />
<col width="10%" />
<col width="5%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">rlpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>group</td>
<td>10</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>rlp</td>
<td>water</td>
<td>vg_1</td>
<td>0.0001</td>
<td>1.0</td>
<td>3.0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>cap</td>
<td>air/water</td>
<td>vg</td>
<td>0.0001</td>
<td>1.0</td>
<td>3.0</td>
<td>3.0</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td>0.05</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>0</td>
<td>10</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p><strong>Example 5</strong>.  Dual porosity air/water problem.  Water (and air, by default) will both have vg rel perm model, calculated using capillary pressure as the independent variable (see Table 2 for details).  Air/water capillary pressure will be vg.  Fracture rel perm and capillary pressure models are specified, as well.</p>
<table border="1" class="docutils">
<colgroup>
<col width="11%" />
<col width="13%" />
<col width="9%" />
<col width="9%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
<col width="9%" />
<col width="3%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">rlpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>group</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>rlp</td>
<td>water</td>
<td>vg_cap</td>
<td>0.0212</td>
<td>1.0</td>
<td>1.62</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>fracture</td>
<td>0.03</td>
<td>1.0</td>
<td>3.00</td>
<td>4.06e-09</td>
<td>2.04e-18</td>
<td>2.93e-04</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>cap</td>
<td>air/water</td>
<td>vg</td>
<td>0.0212</td>
<td>1.0</td>
<td>0.00715</td>
<td>1.62</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td>0.0312</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>fracture</td>
<td>0.03</td>
<td>1.0</td>
<td>12.05</td>
<td>3.00</td>
<td>20.0</td>
<td>0.0001</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>group</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>rlp</td>
<td>water</td>
<td>vg_cap</td>
<td>0.154</td>
<td>1.0</td>
<td>0.371</td>
<td>2.37</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>fracture</td>
<td>0.03</td>
<td>1.0</td>
<td>13.72</td>
<td>3.00</td>
<td>7.14e-09</td>
<td>2.51e-18</td>
<td>9.27e-05</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>cap</td>
<td>air/water</td>
<td>vg_cap</td>
<td>0.154</td>
<td>1.0</td>
<td>0.371</td>
<td>2.37</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td>0.164</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>fracture</td>
<td>0.03</td>
<td>1.0</td>
<td>13.72</td>
<td>3.00</td>
<td>20.0</td>
<td>0.0001</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>-1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>-2</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p><strong>Example 6</strong>.  Water/vapor problem. Group 1 applies a corey relative permeability function (also applied to vapor, by default). cp=0.
Group 2 applies interpolated values from a table, read from an external file “doe_rlpm.table”.  (note:  this file can be generated by FEHM
by using the keyword ‘rel’ in the hist macro)</p>
<table border="1" class="docutils">
<colgroup>
<col width="40%" />
<col width="24%" />
<col width="13%" />
<col width="9%" />
<col width="9%" />
<col width="5%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">rlpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>group</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>rlp</td>
<td>water</td>
<td>corey</td>
<td>0.3</td>
<td>0.1</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>group</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>table</td>
<td>water/vapor</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>file</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>input/doe_rlpm.table</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>140</td>
<td>1</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>For example, the file doe_rlpm.table would contain an arbitrary number header rows (each header row must contain a character in the first column)
followed by an arbitrary number of lines each containing the following information:
saturation, relative permeability (wetting phase), relative permeability (non-wetting phase),and capillary pressure (MPa).</p>
<table border="1" class="docutils">
<colgroup>
<col width="55%" />
<col width="17%" />
<col width="17%" />
<col width="12%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">FEHM V3.00pgi64 10-10-20 QA:NA 10/20/2010 14:23:18</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><strong>* DOE Code Comparison Project, Problem 5, Case A *</strong></td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>Relative permeability and Capillary pressure</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>“Saturation” “Liquid” “Vapor” “Capillary pressure”</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>0.00000000</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.500000000E-01</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.100000000</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.150000000</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.200000000</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.250000000</td>
<td>0.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.300000000</td>
<td>0.732682696E-64</td>
<td>1.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.350000000</td>
<td>0.482253086E-04</td>
<td>0.834442515</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.400000000</td>
<td>0.771604938E-03</td>
<td>0.675154321</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.450000000</td>
<td>0.390625000E-02</td>
<td>0.527343750</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.500000000</td>
<td>0.123456790E-01</td>
<td>0.395061728</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.550000000</td>
<td>0.301408179E-01</td>
<td>0.281201775</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.600000000</td>
<td>0.625000000E-01</td>
<td>0.187500000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.650000000</td>
<td>0.115788966</td>
<td>0.114535108</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.700000000</td>
<td>0.197530864</td>
<td>0.617283951E-01</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.750000000</td>
<td>0.316406250</td>
<td>0.273437500E-01</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.800000000</td>
<td>0.482253086</td>
<td>0.848765432E-02</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.850000000</td>
<td>0.706066744</td>
<td>0.110918210E-02</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>0.900000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-even"><td>0.950000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
<td>0.00000000</td>
</tr>
<tr class="row-odd"><td>1.00000000</td>
<td>1.00000000</td>
<td>0.00000000</td>
<td>0.00000000</td>
</tr>
</tbody>
</table>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
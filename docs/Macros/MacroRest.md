---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rest</span></code><a class="headerlink" href="#rest" title="Permalink to this headline">¶</a></h1>

<p>The <span class="pre">rest</span> macro controls the content and format of the initial condition file.</p> 

<p>FEHM writes a restart file for each run. The restart output file name may be given in the input
control file or as terminal input, or if unspecified will default to fehmn.fin. The restart file is used on a subsequent runs by providing the name of the file.
  
It is recommended that the restart input file name be modified to avoid confusion with the restart output file. For example, by
changing the suffix to .ini, the default restart output file, fehmn.fin would be renamed fehmn.ini. Values from
the restart file will overwrite any variable initialization prescribed in the input file. The initial
time of simulation will also be taken from the restart file unless specified in the macro time
input.
</p>
<p>Keywords can be used to define the restart and final condition files. Use 'rsti:' filename and/or 'rsto:' filename. For example, <pre>rsti:  input/doe.ini</pre> will use the file "input/doe.ini" as a restart file.</p>
<br>

<ul class="simple">
<li>Group 1 -     CHDUM</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="11%" />
<col width="75%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>CHDUM</td>
<td>character*80</td>
<td>Keyword(s).  Keywords are entered one per line and terminated with ‘end’ or a blank line.</td>
</tr>
</tbody>
</table>
<p>Keywords are entered one per line and terminated with ‘end’ or a blank line.</p>
<p>Two keywords ‘read’ and ‘write’ are followed by a list of variables; all others stand alone.</p>
<p>Valid keywords (case insensitive) are:</p>
<ul class="simple">
<li>Keywords that control format:<ul>
<li>‘ascii’ - both read and write restart files are ascii (default)</li>
<li>‘binary’ - both read and write restart files are unformatted</li>
<li>‘rbinary’ - unformatted read restart file</li>
<li>‘wbinary’ - unformatted write restart file</li>
<li>‘old’ use old format (input / output format and content is hard-wired)</li>
<li>‘new’ use new format  (input / output format and content is controlled by restart macro keywords)</li>
</ul>
</li>
<li>Keywords that control flux output:<ul>
<li>‘noflux’ - do not output flux</li>
<li>‘flux’ - output liquid and vapor flux</li>
<li>‘lflux’ -  output liquid flux</li>
<li>‘vflux’ - output vapor flux</li>
</ul>
</li>
<li>Keywords that control the list of variables to read / write
(default is to read all variables in restart file; write
all variables in current simulation)<ul>
<li>‘read’ (followed by list of variables, on same line)</li>
<li>‘write’ (followed by list of variables, on same line)</li>
</ul>
</li>
</ul>
<p>Possible read/write variables:</p>
<ul class="simple">
<li>none</li>
<li>all</li>
<li>temp</li>
<li>pres</li>
<li>poro</li>
<li>trac</li>
<li>ptrk</li>
<li>gasp</li>
<li>pini</li>
<li>saturation</li>
<li>co2</li>
<li>mass</li>
<li>disp (disx, disy, disz)</li>
<li>strs or stre or strs (strx, stry, strz, stxy, stxz, styz)</li>
</ul>
<br>
<p>The following is an example of rest. In this example restart data will be written
to an unformatted file and liquid flux will be output. If a read restart file
is used it will be in ascii format and if liquid flux data is present in the
file it will be read.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>rest</td>
</tr>
<tr class="row-even"><td>wbinary</td>
</tr>
<tr class="row-odd"><td>lflux</td>
</tr>
<tr class="row-even"><td>&#160;</td>
</tr>
</tbody>
</table>
<br>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>

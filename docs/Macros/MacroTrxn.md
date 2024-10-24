---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">trxn</span></code><a class="headerlink" href="#trxn" title="Permalink to this headline">¶</a></h1>
<div class="section" id="overview">
<h2>Overview<a class="headerlink" href="#overview" title="Permalink to this headline">¶</a></h2>
<p>The <code class="docutils literal notranslate"><span class="pre">trxn</span></code> macro is designed to replace the <code class="docutils literal notranslate"><span class="pre">trac</span></code> (and optionally <code class="docutils literal notranslate"><span class="pre">rxn</span></code>) macros with a more convenient, user-friendly input format.  Rather than reading input as sets of nameless numbers, the <code class="docutils literal notranslate"><span class="pre">trxn</span></code> macro reads several “blocks” of data in which parameters are defined.  The blocks may be specified in any order (with the exception of the <code class="docutils literal notranslate"><span class="pre">ctrl</span></code> and <code class="docutils literal notranslate"><span class="pre">lookup</span></code> blocks), and all blocks are optional, except for the <code class="docutils literal notranslate"><span class="pre">header</span></code> and <code class="docutils literal notranslate"><span class="pre">comp</span></code> blocks.  Any parameters that are not specified will be given default values (usually zero).</p>
<p>The <code class="docutils literal notranslate"><span class="pre">trxn</span></code> macro relies heavily on zones, and uses zones for applying all variables that can vary by node.  For this reason, a <code class="docutils literal notranslate"><span class="pre">zone</span></code> macro must be supplied in the input file before <cite>trxn</cite> is read.  A <code class="docutils literal notranslate"><span class="pre">time</span></code> macro must also be given before <code class="docutils literal notranslate"><span class="pre">trxn</span></code>.</p>
<p>If the pound character (“#”) appears on a line, everything after it on that line will be treated as a comment.  Lines in which the first character is a pound sign are treated as blank lines.  Blocks must not contain any blank or commented lines (although they can contain comments that do not start at the beginning of the line), and blocks must be separated by at least one blank or commented line.  In the <code class="docutils literal notranslate"><span class="pre">comp</span></code>, <code class="docutils literal notranslate"><span class="pre">water</span></code>, <code class="docutils literal notranslate"><span class="pre">rock</span></code>, <code class="docutils literal notranslate"><span class="pre">gas</span></code>, <code class="docutils literal notranslate"><span class="pre">disp</span></code>, <code class="docutils literal notranslate"><span class="pre">sorp</span></code>, and <code class="docutils literal notranslate"><span class="pre">assign</span></code> blocks, entire columns can be commented out.  If an asterisk (“*”) (separated from surrounding tokens by whitespace) appears in the column header line of one of these blocks, the contents of every column to the right of the asterisk will be ignored.  An entire block (until the next blank line) can be skipped by placing the keyword <code class="docutils literal notranslate"><span class="pre">null</span></code> before the block name.</p>
</div>
<div class="section" id="trxn-blocks">
<h2>trxn blocks<a class="headerlink" href="#trxn-blocks" title="Permalink to this headline">¶</a></h2>
<p>The blocks are as follows:</p>
</div>
<div class="section" id="ctrl">
<h2><code class="docutils literal notranslate"><span class="pre">ctrl</span></code><a class="headerlink" href="#ctrl" title="Permalink to this headline">¶</a></h2>
<p>This block contains control parameters, all on the same line.  If it is supplied, it must be the first block in trxn.  Its options are as follows:</p>
<table border="1" class="docutils">
<colgroup>
<col width="9%" />
<col width="91%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Option</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">rxnon</span></code></td>
<td>Enables reactions, which are by default disabled.  Reactions will not occur in a simulation if <cite>rxnon</cite> is not given, even if reaction-related blocks are specified.</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">co2_couple</span></code></td>
<td>Enables CO2 coupling for <cite>rxn</cite>.</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">debug</span></code></td>
<td>Enables the output of debugging information.</td>
</tr>
</tbody>
</table>
<p>Below is an example of the <cite>ctrl</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">ctrl</span> <span class="n">co2_couple</span> <span class="n">rxnon</span>
</pre></div>
</div>
<p>In this example, CO2 coupling and reactions are enabled.</p>
</div>
<div class="section" id="include">
<h2><code class="docutils literal notranslate"><span class="pre">include</span></code><a class="headerlink" href="#include" title="Permalink to this headline">¶</a></h2>
<p>The <cite>include</cite> block is designed to facilitate the construction of “libraries” of rock/water/gas types, component properties, dispersion models, etc. that can be shared between simulations.  It allows blocks to be included from external files as though these files’ contents were placed at the location of the <cite>include</cite> statement.  External include libraries have the exact same syntax as the standard <cite>trxn</cite> macro.  The first non-blank (and non-commented) line of the include file should be <cite>trxn library</cite>.  Following this, as many blocks as desired may be supplied.
The <cite>end trxn</cite> line should be given at the end of the library file;
otherwise, reading will terminate early.
Multiple <cite>include</cite> statements referencing several different files may be used in
the same <cite>trxn</cite> macro.  Libraries may include other libraries.</p>
<p>The syntax is as follows:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">include</span> <span class="o">/</span><span class="n">scratch</span><span class="o">/</span><span class="n">nts</span><span class="o">/</span><span class="n">ms</span><span class="o">/</span><span class="n">trxn</span><span class="o">/</span><span class="n">test</span><span class="o">-</span><span class="n">problems</span><span class="o">/</span><span class="n">standard</span><span class="o">.</span><span class="n">trxl</span>
</pre></div>
</div>
<p>This will include the library located at <code class="docutils literal notranslate"><span class="pre">/scratch/nts/ms/trxn/test-problems/standard.trxl</span></code> at the current location in the input file.  As an example, a simple library might look something like this:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># Standard library for use with all test problems</span>

<span class="n">trxn</span> <span class="n">library</span>

<span class="n">header</span>
<span class="mi">0</span>    <span class="mi">1</span>       <span class="mf">1e-6</span>    <span class="mi">1</span>
<span class="mi">0</span>    <span class="mi">100000</span>  <span class="mi">100000</span>  <span class="mi">100000</span>
<span class="mi">5</span>    <span class="mi">2</span>       <span class="mi">1</span>       <span class="mi">1</span>       <span class="mi">0</span>
<span class="n">iskip</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">rsdmax</span><span class="o">=</span><span class="mf">1e-9</span>

<span class="n">disp</span> <span class="n">lx</span>      <span class="n">ly</span>      <span class="n">lz</span>      <span class="n">vx</span>      <span class="n">vy</span>      <span class="n">vz</span>
<span class="n">std_x</span>        <span class="mf">0.1</span>     <span class="mf">1e-9</span>    <span class="mf">1e-9</span>    <span class="mf">0.1</span>     <span class="mf">1e-9</span>    <span class="mf">1e-9</span>
<span class="n">std_y</span>        <span class="mf">1e-9</span>    <span class="mf">0.1</span>     <span class="mf">1e-9</span>    <span class="mf">1e-9</span>    <span class="mf">0.1</span>     <span class="mf">1e-9</span>
<span class="n">std_z</span>        <span class="mf">1e-9</span>    <span class="mf">1e-9</span>    <span class="mf">0.1</span>     <span class="mf">1e-9</span>    <span class="mf">1e-9</span>    <span class="mf">0.1</span>

<span class="n">sorp</span> <span class="n">ltype</span>   <span class="n">a1l</span>     <span class="n">a2l</span>     <span class="n">bl</span>      <span class="n">vtype</span>   <span class="n">a1v</span>     <span class="n">a2v</span>     <span class="n">bv</span>
<span class="o">.</span><span class="n">std</span>
<span class="o">*</span>    <span class="n">con</span>     <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>       <span class="n">con</span>     <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>

<span class="n">diff</span> <span class="n">l</span><span class="o">=</span><span class="mf">1e-9</span><span class="p">,</span> <span class="n">v</span><span class="o">=</span><span class="mf">1e-9</span>

<span class="n">end</span> <span class="n">trxn</span>
</pre></div>
</div>
<p>All problems that include this library will be given the standard values for <cite>header</cite> and liquid and vapor diffusion coefficients of 1×10^-9^, and will have one sorption model (“std”) and three dispersion models (“std_x”, “std_y”, and “std_z”) available.</p>
<p>Caution should be used when including libraries to prevent multiple definitions of any blocks, as this may cause unexpected behavior.  <cite>trxn</cite> will print a warning message to the error file if a block is supplied more than once, and these warnings should be considered.  Also, care should be taken to avoid placing any problem-specific parameters in the libraries.</p>
</div>
<div class="section" id="header">
<h2><code class="docutils literal notranslate"><span class="pre">header</span></code><a class="headerlink" href="#header" title="Permalink to this headline">¶</a></h2>
<p>This block contains basic constant values, and is the first three lines of the <cite>trac</cite> macro copied verbatim.</p>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="87%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Line</th>
<th class="head">Variables</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td><code class="docutils literal notranslate"><span class="pre">AN0</span> <span class="pre">AWC</span> <span class="pre">EPC</span> <span class="pre">UPWGTA</span></code></td>
</tr>
<tr class="row-odd"><td>2</td>
<td><code class="docutils literal notranslate"><span class="pre">DAYCS</span> <span class="pre">DAYCF</span> <span class="pre">DAYHF</span> <span class="pre">DAYHS</span></code></td>
</tr>
<tr class="row-even"><td>3</td>
<td><code class="docutils literal notranslate"><span class="pre">IACCMX</span> <span class="pre">DAYCM</span> <span class="pre">DAYCMM</span> <span class="pre">DAYCMX</span> <span class="pre">NPRTTRC</span></code></td>
</tr>
<tr class="row-odd"><td>4</td>
<td>Optional parameters</td>
</tr>
</tbody>
</table>
<p>Please refer to the <cite>trac</cite> section of the FEHM User’s Manual for the meanings of the parameters on the first three lines.  On the fourth line, optional parameters may be defined, using the form <cite>variable=value</cite> (no spaces), with commas between these pairs.  The variables that may be set here are <cite>ISKIP</cite>, <cite>RSDMAX</cite>, and <cite>STRAC_MAX</cite>.  Please refer to the <cite>rxn</cite> section of the FEHM User’s Manual for the meanings of these variables.  (<cite>STRAC_MAX</cite> is the maximum allowable saturation for vapor and Henry’s Law species.  When using <cite>trac</cite>, it is read if provided in the input file; however, it is not documented in the <cite>trac</cite> portion of the FEHM User’s Manual.)  <cite>ISKIP</cite> and <cite>RSDMAX</cite> are used only if reactions are enabled (see <cite>rxn</cite> below).  If these optional variables are omitted, they are given the values shown in the example below.</p>
<p>Below is an example of the <cite>header</cite> block.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">header</span>
<span class="mf">0.0</span>  <span class="mf">1.0</span>     <span class="mf">1.0e-6</span>  <span class="mf">1.0</span>                     <span class="c1"># ANO, AWC, EPC, UPWGTA</span>
<span class="mf">0.0</span>  <span class="mf">1e20</span>    <span class="mf">1e20</span>    <span class="mf">1e20</span>                    <span class="c1"># DAYCS, DAYCF, DAYHF, DAYHS</span>
<span class="mi">10</span>   <span class="mf">2.0</span>     <span class="mf">1.0</span>     <span class="mf">150000.0</span>        <span class="mi">1</span>       <span class="c1"># IACCMX, DAYCM, DAYCMM, DAYCMX, NPRTTRC</span>
<span class="n">iskip</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">rsdmax</span><span class="o">=</span><span class="mf">1e-9</span><span class="p">,</span> <span class="n">strac_max</span><span class="o">=</span><span class="mf">0.99</span>
</pre></div>
</div>
<p>Please note that the keyword <code class="docutils literal notranslate"><span class="pre">userc</span></code> is not supported in the <code class="docutils literal notranslate"><span class="pre">header</span></code> block.  For <code class="docutils literal notranslate"><span class="pre">userc</span></code> support, please refer to the <cite>userc</cite> block below.  Also, unlike in <code class="docutils literal notranslate"><span class="pre">trac</span></code>, <code class="docutils literal notranslate"><span class="pre">NPRTTRC</span></code> may not be omitted from the header.</p>
</div>
<div class="section" id="userc">
<h2><code class="docutils literal notranslate"><span class="pre">userc</span></code><a class="headerlink" href="#userc" title="Permalink to this headline">¶</a></h2>
<p>This block invokes the solute transport user subroutine as specified in the <code class="docutils literal notranslate"><span class="pre">trac</span></code> section of the user’s manual.  On the same line as the macro name should be given the path to a file containing userc parameters.  See the <cite>trac</cite> section of the user’s manual for more information on the <cite>userc</cite> input format.</p>
<p>Below is an example of <code class="docutils literal notranslate"><span class="pre">userc</span></code>:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">userc</span> <span class="nb">input</span><span class="o">/</span><span class="n">userc</span><span class="o">.</span><span class="n">dat</span>
</pre></div>
</div>
<p>In this example, the <code class="docutils literal notranslate"><span class="pre">userc</span></code> subroutine is called with an input file located at <code class="docutils literal notranslate"><span class="pre">input/userc.dat</span></code>.</p>
</div>
<div class="section" id="comp">
<h2><code class="docutils literal notranslate"><span class="pre">comp</span></code><a class="headerlink" href="#comp" title="Permalink to this headline">¶</a></h2>
<p>The <code class="docutils literal notranslate"><span class="pre">comp</span></code> block is used to define each component present in the simulation.  (A component is any group of compounds, ions, etc., all of the same phase, that the user wishes to be treated as a single entity by the tracer solver.)  It contains one line for each component, and each line consists of a phase designation for that component and the name of the component.  The phase designation is one of “aqueous”, “solid”, “gas”, and “henry”; these may be shortened to the first character to save time.</p>
<p>If kinetic reactions (<code class="docutils literal notranslate"><span class="pre">rxn</span></code> macro) are being simulated, two additional columns may be included.  The <cite>master</cite> column indicates the “master species” for each aqueous or Henry’s Law component.  These master species are arbitrarily chosen forms of the components, by convention the form that is expected to dominate in reactions.  The third column, <cite>guess</cite>, allows the user to specify a guess for the initial uncomplexed concentration of each aqueous and Henry’s Law component.  This is not necessary unless the chemical speciation solver has difficulty converging with the default value of 1×10^-9^, in which case specifying a more representative value may help.  If this is not necessary, leave the column out entirely, or place asterisks in the rows of components that do not need help converging.</p>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">comp</span></code> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">comp</span>         <span class="n">master</span>  <span class="n">guess</span>
<span class="n">aqueous</span> <span class="n">H</span>    <span class="n">H</span><span class="o">+</span>      <span class="o">*</span>
<span class="n">aq</span> <span class="n">C_a</span>               <span class="n">HCO3</span><span class="o">-</span>   <span class="mf">8.6e-6</span>
<span class="n">a</span> <span class="n">Na_a</span>               <span class="n">Na</span><span class="o">+</span>     <span class="o">*</span>
<span class="n">a</span> <span class="n">Ca_a</span>               <span class="n">Ca</span><span class="o">++</span>    <span class="o">*</span>
<span class="n">a</span> <span class="n">C_a2</span>               <span class="n">CO3</span><span class="o">--</span>   <span class="o">*</span>
<span class="n">a</span> <span class="n">Cl_a</span>               <span class="n">Cl</span><span class="o">-</span>     <span class="mf">2.0e-9</span>
<span class="n">a</span> <span class="n">U238</span>               <span class="n">U02</span>     <span class="o">*</span>
<span class="n">a</span> <span class="n">Th234</span>              <span class="n">ThO2</span>    <span class="o">*</span>
<span class="n">solid</span> <span class="n">AlO3</span>   <span class="o">*</span>       <span class="o">*</span>
<span class="n">s</span> <span class="n">NaCl</span>               <span class="o">*</span>       <span class="o">*</span>
<span class="n">s</span> <span class="n">NaHCO3</span>     <span class="o">*</span>       <span class="o">*</span>
<span class="n">s</span> <span class="n">CaCl2</span>              <span class="o">*</span>       <span class="o">*</span>
<span class="n">gas</span> <span class="n">O2</span>               <span class="o">*</span>       <span class="o">*</span>
<span class="n">g</span> <span class="n">N2</span>         <span class="o">*</span>       <span class="o">*</span>
<span class="n">henry</span> <span class="n">C_h2</span>   <span class="n">C6H6</span>    <span class="mf">1.2e-8</span>
<span class="n">h</span> <span class="n">Cl_h</span>               <span class="n">Cl2</span>     <span class="o">*</span>
<span class="n">h</span> <span class="n">C_h</span>                <span class="n">CO2</span>     <span class="o">*</span>
</pre></div>
</div>
<p>In this example, there are 17 components.  “H”, “C_a”, “Na_a”, “Ca_a”, “C_a2”,
“Cl_a”, “U238”, and “Th234” are aqueous; “AlO3”, “NaCl”, “NaHCO3”, and “CaCl2”
are solid; “O2” and “N2” are gaseous, and “C_h”, “Cl_h”, and “C_h2” may be liquid
or gas according to Henry’s Law.</p>
</div>
<div class="section" id="water-rock-gas">
<h2><code class="docutils literal notranslate"><span class="pre">water</span></code>, <code class="docutils literal notranslate"><span class="pre">rock</span></code>, <code class="docutils literal notranslate"><span class="pre">gas</span></code><a class="headerlink" href="#water-rock-gas" title="Permalink to this headline">¶</a></h2>
<dl class="docutils">
<dt>These blocks are identical in form, and they are used to assign concentrations in</dt>
<dd>the simulation for components of different states.  These blocks allow the user to
specify different “water types”, “rock types”, and “gas types”, which may consist
of different combinations of components specified in <cite>comp</cite> in different concentrations.
See the <cite>moles</cite> block below for an alternative input format.</dd>
<dt>On the same line as the block name, the names of components specified in <cite>comp</cite></dt>
<dd>are placed, separated by tabs or spaces.  These are column headers.  Below this line,
one line is given to each “type” desired.  Each line consists of the name of the type,
followed by numbers representing the concentrations of each of the components given
in the column headers in that type.  If the columns are separated by tabs,
this layout forms a neat table with types down the left side, components across
the top, and the concentration of each component in each type at the intersections
of the rows and columns.</dd>
</dl>
<p>Only aqueous and Henry’s Law components may be included in the <cite>water</cite> block,
only solid components in the <cite>rock</cite> block, and only gaseous and Henry’s Law
components in the <cite>gas</cite> block.  In the <cite>water</cite> block, a special column header,
<cite>pH</cite>, may be included.  This is the same as heading the column with “H”
(and may only be done if “H” is specified in <cite>comp</cite> and is aqueous), but
allows the user to enter H+ concentration in terms of pH rather than molarity.
If a concentration in the grid is negative, it is assumed that the value entered
is the base-ten logarithm of the actual concentration, and the value is adjusted
accordingly.</p>
<p>An asterisk (“*”) in a space where a number is expected is the same as a 0. If a
component is omitted entirely from the table, it is assumed that that
component is not present in the simulation.
The unit for all numbers in these tables is molal, with the exception of the
<cite>pH</cite> column, if present.</p>
<p>A negative value in the grid for a water or gas inflow type indicates that
the concentration of that solute will be held constant at inflow nodes of that water type.
Negative values in initial condition types should be avoided.</p>
<p>An example of each of the four block types, consistent with the sample <cite>comp</cite>
block above, is given below:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">water</span>        <span class="n">pH</span>      <span class="n">C_a</span>     <span class="n">Na_a</span>    <span class="n">C_a2</span>    <span class="n">C_h2</span>
<span class="n">wt1</span>  <span class="mf">7.43</span>    <span class="mf">1e-6</span>    <span class="mf">1.34e-5</span> <span class="mf">1e-10</span>   <span class="mf">0.3</span>
<span class="n">wt2</span>  <span class="mf">5.4</span>     <span class="mf">1e-3</span>    <span class="mf">0.002</span>   <span class="mf">2.4e-6</span>  <span class="mf">2.3</span>
<span class="n">wt3</span>  <span class="mf">2.3</span>     <span class="mi">0</span>       <span class="mi">10</span>      <span class="mi">0</span>       <span class="mi">0</span>

<span class="n">rock</span> <span class="n">AlO3</span>    <span class="n">NaCl</span>    <span class="n">NaHCO3</span>
<span class="n">tuff</span> <span class="mf">3.4e-2</span>  <span class="mf">1.24</span>    <span class="mf">1.6e-6</span>
<span class="n">granite</span>      <span class="mf">3e-5</span>    <span class="mf">1.5e-2</span>  <span class="mf">6e-5</span>
<span class="n">clay</span> <span class="mf">0.3</span>     <span class="mf">3.9e-3</span>  <span class="mf">0.03</span>

<span class="n">gas</span>  <span class="n">O2</span>      <span class="n">N2</span>      <span class="n">Cl_h</span>    <span class="n">C_h</span>
<span class="n">vt1</span>  <span class="mf">0.23</span>    <span class="mf">0.10</span>    <span class="mf">1.24</span>    <span class="mi">2</span>
<span class="n">vt2</span>  <span class="mf">1.02</span>    <span class="mf">0.012</span>   <span class="mf">0.2</span>     <span class="mi">0</span>
</pre></div>
</div>
<p>In this example, there are three water types (“wt1”, “wt2”, and “wt3”), three rock
types (“tuff”, “granite”, and “clay”), and two vapor types (“vt1” and “vt2”).
Water type “wt1” has a pH of 7.3, an HCO3- concentration of 1×10^-6^ //m//, an Na+
concentration of 1.34×10^-5^ //m//, a CO3– concentration of 1×10^-10^ //m//,
and a C6H6 concentration of 0.3 //m//.</p>
</div>
<div class="section" id="print">
<h2><code class="docutils literal notranslate"><span class="pre">print</span></code><a class="headerlink" href="#print" title="Permalink to this headline">¶</a></h2>
<p>The <cite>print</cite> block allows the user to specify for which components and complexes information is to be printed to tracer output(<cite>.trc</cite>) files.  This block occupies only one line.  After the keyword <cite>print</cite> is given a list of aqueous component and complex names that are to be printed, delimited by spaces, tabs, or commas.  The keywords <cite>all</cite> and <cite>none</cite> may be given instead of the list, specifying that information is to be printed for all or none of the aqueous components and complexes, respectively.  The default action, if no <cite>print</cite> block is specified, is the same as specifying <cite>print all</cite>.  This information is only used if reactions are enabled.</p>
<p>Below is an example of the <cite>print</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="nb">print</span> <span class="n">UO2</span> <span class="n">HCO3</span><span class="o">-</span> <span class="n">Na_a</span> <span class="n">Ca_a</span>
</pre></div>
</div>
<p>In this example, information will be printed only about UO2, HCO3-, Na_a, and Ca_a.</p>
</div>
<div class="section" id="moles">
<h2><code class="docutils literal notranslate"><span class="pre">moles</span></code><a class="headerlink" href="#moles" title="Permalink to this headline">¶</a></h2>
<p>This block allows the user to specify initial solute conditions as the total number of moles contained in each zone.  FEHM will distribute the specified number of moles evenly throughout the volume of the zone.  The <cite>moles</cite> block conflicts with the <cite>water</cite>, <cite>rock</cite>, and <cite>gas</cite> columns in the <cite>assign</cite> block, and they cannot both appear in the same <cite>trxn</cite> macro.</p>
<p>The <cite>moles</cite> block allows the use of custom zone definitions.  Unlike standard zones, these zones are permitted to overlap, and the concentrations at nodes in overlapping zones become cumulative.  If these custom zones are to be used, a separate file must be created, containing a valid zone macro that defines the desired zones.  Due to limited functionality in the <cite>trac</cite> macro, all zones in this macro must be specified using the <cite>nnum</cite> method and zones must be numbered rather than named.  (See the manual section on the <cite>zone</cite> macro for more information.)  Alternatively, the zones defined in previous <cite>zone</cite> macros may be used.</p>
<p>On the same line as the <cite>moles</cite> keyword in this block, the path to the alternate zone file should be provided.  If the standard zone definitions are to be used, leave the path blank or explicitly turn the alternate zone processing off by providing an asterisk in the place of the path.  The next lines form a table with component names from <cite>comp</cite> across the top and zone numbers down the left-hand side.  At the intersections of rows and columns is placed the number of moles of that component initially in that zone.  Any omitted zones or components are assigned a zero concentration.</p>
<p>Below is an example of the <cite>moles</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">moles</span> <span class="nb">input</span><span class="o">/</span><span class="n">moles</span><span class="o">.</span><span class="n">zone</span>
     <span class="n">C_a</span>     <span class="n">Na_a</span>    <span class="n">C_a2</span>    <span class="n">C_h2</span>
<span class="mi">2</span>    <span class="mi">0</span>       <span class="mf">1.6</span>     <span class="mf">0.02</span>    <span class="mi">1</span>
<span class="mi">3</span>    <span class="mi">0</span>       <span class="mf">1.7</span>     <span class="mf">0.6</span>     <span class="mf">0.89</span>
<span class="mi">5</span>    <span class="mi">0</span>       <span class="mi">2</span>       <span class="mf">0.045</span>   <span class="mf">2.2</span>
<span class="mi">4</span>    <span class="mf">1.2</span>     <span class="mi">0</span>       <span class="mf">0.102</span>   <span class="mf">1.61</span>
</pre></div>
</div>
<p>In this example, zone definitions are loaded from the file <cite>input/moles.zone</cite>, which remain in effect for the rest of the <cite>moles</cite> macro.  Zone 2 contains 1.6 moles of Na_a, 0.02 moles of C_a2, and 1 mole of C_h2 distributed evenly throughout it.</p>
</div>
<div class="section" id="hparam">
<h2><code class="docutils literal notranslate"><span class="pre">hparam</span></code><a class="headerlink" href="#hparam" title="Permalink to this headline">¶</a></h2>
<p>The <code class="docutils literal notranslate"><span class="pre">hparam</span></code> block is used to assign Henry’s Law parameters to different Henry’s Law components.  Each Henry’s Law component appearing in <cite>comp</cite> must be included in this block.  Below the block title <cite>hparam</cite>, each Henry’s Law component is given one line.  Each of these lines consists of the name of the Henry’s Law component, the model that will be used to simulate that component, and then several <cite>parameter=value</cite> pairs (separated by commas), specifying parameters for each model.  The models and their parameters are as follows:</p>
<table border="1" class="docutils">
<colgroup>
<col width="12%" />
<col width="88%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Option</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><cite>1</cite> or <cite>hoff</cite></td>
<td>van’t Hoff model.  Parameters:  <span class="math notranslate nohighlight">\(ah - A_H; dhh - \delta H_H\)</span></td>
</tr>
<tr class="row-odd"><td><cite>2</cite> or <cite>multi</cite></td>
<td>Multi-parameter fit.  Parameters: <span class="math notranslate nohighlight">\(ahn - A_Hn\)</span> for each n from 1 to 5</td>
</tr>
<tr class="row-even"><td><cite>3</cite> or <cite>wvp</cite></td>
<td>Use water vapor pressure.  Parameters: <span class="math notranslate nohighlight">\(ah-A_H\)</span>; hh - Henry’s constant modifier <span class="math notranslate nohighlight">\(H=P_{wv} \cdot \Delta H_H\)</span></td>
</tr>
</tbody>
</table>
<p>Below is an example of an <cite>hparam</cite> block that uses all three of the above models:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">hparam</span>
<span class="n">C6H6</span> <span class="n">hoff</span>    <span class="n">ah</span><span class="o">=</span><span class="mf">0.34</span><span class="p">,</span> <span class="n">dhh</span><span class="o">=</span><span class="mf">3.2</span>
<span class="n">Cl2</span>  <span class="n">multi</span>   <span class="n">ah1</span><span class="o">=</span><span class="mf">4.12</span><span class="p">,</span> <span class="n">ah2</span><span class="o">=</span><span class="mf">2.4</span><span class="p">,</span> <span class="n">ah3</span><span class="o">=</span><span class="mi">4</span><span class="p">,</span> <span class="n">ah4</span><span class="o">=</span><span class="mf">0.18</span><span class="p">,</span> <span class="n">ah5</span><span class="o">=</span><span class="mf">1.1</span>
<span class="n">CO2</span>  <span class="n">wvp</span>     <span class="n">hh</span><span class="o">=</span><span class="mf">2.8</span><span class="p">,</span> <span class="n">ah</span><span class="o">=</span><span class="mf">0.04</span>
</pre></div>
</div>
<p>In this example, Henry’s Law component C6H6 is simulated using the van’t Hoff model.  Its A,,H,, value is 0.34, and its ΔH,,H,, value is 3.2.</p>
</div>
<div class="section" id="diff">
<h2><code class="docutils literal notranslate"><span class="pre">diff</span></code><a class="headerlink" href="#diff" title="Permalink to this headline">¶</a></h2>
<p>This block sets the molecular diffusion coefficients and diffusion models for liquid
and vapor components.  Due to technical restrictions imposed by FEHM’s solver,
there can only be one liquid and one vapor diffusion coefficient for the entire model.
The <code class="docutils literal notranslate"><span class="pre">diff</span></code> block consists of one line.  After the keyword are up to four name/value pairs, where the name is <cite>l</cite> to set the liquid diffusion coefficient, <cite>v</cite> to set the vapor diffusion coefficient, <cite>lm</cite> to set the liquid diffusion model, or <cite>vm</cite> to set the liquid diffusion model.  The possible values for the diffusion models are as follows:</p>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="88%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Option</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><cite>0</cite> or <cite>con</cite></td>
<td>Constant molecular diffusion coefficient</td>
</tr>
<tr class="row-odd"><td><cite>1</cite> or <cite>mq</cite></td>
<td>Millington-Quirk diffusion model</td>
</tr>
<tr class="row-even"><td><cite>2</cite> or <cite>cw</cite></td>
<td>Conca-Wright diffusion model for liquid, or alternate Millington-Quirk model for vapor</td>
</tr>
<tr class="row-odd"><td><cite>3</cite> or <cite>adif</cite></td>
<td>Diffusion model calculated as a function of pressure and temperature using tortuosity from <cite>adif</cite> macro</td>
</tr>
</tbody>
</table>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">diff</span></code> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">diff</span> <span class="n">l</span><span class="o">=</span><span class="mf">1e-9</span><span class="p">,</span> <span class="n">v</span><span class="o">=</span><span class="mf">1.6e-8</span>        <span class="n">lm</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">vm</span><span class="o">=</span><span class="n">con</span>
</pre></div>
</div>
<p>In this example, the liquid diffusion coefficient is 1×10^-9^, and the vapor diffusion coefficient is 1.6×10^-8^.  Both liquid and vapor use a constant molecular diffusion coefficient.</p>
</div>
<div class="section" id="disp">
<h2><code class="docutils literal notranslate"><span class="pre">disp</span></code><a class="headerlink" href="#disp" title="Permalink to this headline">¶</a></h2>
<p>This block is used to set the dispersivity constants of regions of the simulation.  Dispersivity parameters are applied to dispersivity models, which are applied to zones in the <cite>assign</cite> block.  Model names run down the left side of the block, and parameters run across the top.  Dispersivity can be supplied either for the X, Y, and Z directions, or for the longitudinal and transverse directions.  Parameter names consist of two characters:  the first <cite>l</cite> or <cite>v</cite> for liquid or vapor, and the second <cite>x</cite>, <cite>y</cite>, <cite>z</cite>, <cite>l</cite>, or <cite>t</cite> for X, Y, Z, longitudinal, or transverse, respectively.  Parameter names from the two modes of setting dispersion cannot be mixed.</p>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">disp</span></code> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">disp</span> <span class="n">lx</span>      <span class="n">ly</span>      <span class="n">lz</span>      <span class="n">vx</span>      <span class="n">vy</span>      <span class="n">vz</span>
<span class="n">model1</span>       <span class="mf">0.2</span>     <span class="mi">3</span>       <span class="mf">1.5</span>     <span class="mf">0.28</span>    <span class="mf">3.6</span>     <span class="mf">1.92</span>
<span class="n">model2</span>       <span class="mf">0.18</span>    <span class="mf">0.9</span>     <span class="mf">2.361</span>   <span class="mf">1.22</span>    <span class="mf">0.56</span>    <span class="mf">0.58</span>
</pre></div>
</div>
<p>In this example, X/Y/Z dispersivity is set.  There are two models, named <cite>model1</cite> and <cite>model2</cite>.  The liquid dispersivity for <cite>model1</cite> is 0.2 in the X direction, 3 in the Y direction, and 1.5 in the Z direction.  The vapor dispersivity for <cite>model1</cite> is 0.28 in the X direction, 3.6 in the Y direction, and 1.92 in the Z direction.</p>
</div>
<div class="section" id="sorp">
<h2><code class="docutils literal notranslate"><span class="pre">sorp</span></code><a class="headerlink" href="#sorp" title="Permalink to this headline">¶</a></h2>
<p>This block sets adsorption parameters for selected components.  On the same line as the block title should appear any or all of the following column headers:</p>
<table border="1" class="docutils">
<colgroup>
<col width="6%" />
<col width="94%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Header</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><cite>ltype</cite></td>
<td>The type of adsorption model to be used to simulate liquid adsorption.  See below for the possible adsorption models from which to choose.</td>
</tr>
<tr class="row-odd"><td><cite>a1l</cite></td>
<td><span class="math notranslate nohighlight">\(\alpha_1\)</span> parameter in liquid adsorption model.</td>
</tr>
<tr class="row-even"><td><cite>a2l</cite></td>
<td><span class="math notranslate nohighlight">\(\alpha_2\)</span> parameter in liquid adsorption model.</td>
</tr>
<tr class="row-odd"><td><cite>bl</cite></td>
<td><span class="math notranslate nohighlight">\(\beta\)</span> parameter in liquid adsorption model.</td>
</tr>
<tr class="row-even"><td><cite>vtype</cite></td>
<td>The type of adsorption model to be used to simulate vapor adsorption.  See below for the possible adsorption models from which to choose.</td>
</tr>
<tr class="row-odd"><td><cite>a1v</cite></td>
<td><span class="math notranslate nohighlight">\(\alpha_1\)</span> parameter in vapor adsorption model.</td>
</tr>
<tr class="row-even"><td><cite>a2v</cite></td>
<td><span class="math notranslate nohighlight">\(\alpha_2\)</span> parameter in vapor adsorption model.</td>
</tr>
<tr class="row-odd"><td><cite>bv</cite></td>
<td><span class="math notranslate nohighlight">\(\beta\)</span> parameter in vapor adsorption model.</td>
</tr>
</tbody>
</table>
<p>If a column header is not specified, it is assumed to be zero for every component.</p>
<p>The following are the available adsorption models to use for <code class="docutils literal notranslate"><span class="pre">ltype</span></code> and <code class="docutils literal notranslate"><span class="pre">vtype</span></code>:</p>
<table border="1" class="docutils">
<colgroup>
<col width="29%" />
<col width="71%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Option</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><cite>0</cite> or <cite>con</cite></td>
<td>Conservative solute</td>
</tr>
<tr class="row-odd"><td><cite>1</cite> or <cite>lin</cite></td>
<td>Linear sorption isotherm</td>
</tr>
<tr class="row-even"><td><cite>2</cite> or <cite>freu</cite></td>
<td>Freundlich sorption isotherm</td>
</tr>
<tr class="row-odd"><td><cite>3</cite> or <cite>mfreu</cite></td>
<td>Modified Freundlich sorption isotherm</td>
</tr>
<tr class="row-even"><td><cite>4</cite> or <cite>lang</cite></td>
<td>Langmuir sorption isotherm</td>
</tr>
</tbody>
</table>
<p>The <span class="math notranslate nohighlight">\(\alpha_1, \alpha_2,\)</span> and <span class="math notranslate nohighlight">\(\beta\)</span> parameters are used differently
according to the adsorption model chosen:</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="38%" />
<col width="21%" />
<col width="12%" />
<col width="15%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Model</th>
<th class="head">Expression</th>
<th class="head"><span class="math notranslate nohighlight">\(\alpha_1\)</span></th>
<th class="head"><span class="math notranslate nohighlight">\(\alpha_2\)</span></th>
<th class="head"><span class="math notranslate nohighlight">\(\beta\)</span></th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>Linear</td>
<td><span class="math notranslate nohighlight">\(C_r = K_d \cdot C_l\)</span></td>
<td><span class="math notranslate nohighlight">\(K_d\)</span></td>
<td>0</td>
<td>1</td>
</tr>
<tr class="row-odd"><td>Freundlich</td>
<td><span class="math notranslate nohighlight">\(C_r = \Lambda \cdot C_1^\beta\)</span></td>
<td><span class="math notranslate nohighlight">\(\Lambda\)</span></td>
<td>0</td>
<td><span class="math notranslate nohighlight">\(0 &lt; \beta &lt; 1\)</span></td>
</tr>
<tr class="row-even"><td>Modified Freundlich</td>
<td><span class="math notranslate nohighlight">\(C_r / (C_{r,max} - C_r) = \Lambda \cdot C_l^\beta\)</span></td>
<td><span class="math notranslate nohighlight">\(\Lambda \cdot C_{r,max}\)</span></td>
<td><span class="math notranslate nohighlight">\(\Lambda\)</span></td>
<td><span class="math notranslate nohighlight">\(0 &lt; \beta &lt; 1\)</span></td>
</tr>
<tr class="row-odd"><td>Langmuir</td>
<td><span class="math notranslate nohighlight">\(C_r = (r_b \cdot C_l) / (1 + r \cdot C_l)\)</span></td>
<td><span class="math notranslate nohighlight">\(r_b\)</span></td>
<td>r</td>
<td>1</td>
</tr>
</tbody>
</table>
<p>For a more in-depth description of the models used for sorption, please refer to the FEHM Models and Methods Summary.</p>
<p>Models are designated by a line starting with a period (“.”), immediately followed by the name of the model and nothing else.  On the next several lines, components specified in <cite>comp</cite> may be given, with their corresponding parameters for the current model in the correct column.  The same components, in the same order, must be given in each model.  If a component is liquid- or vapor-only, then asterisks should be placed in the columns that do not apply to that component.  An asterisk can also be placed in the column containing component names to indicate that all components that have not yet been explicitly given sorption parameters should use the model on the line with the asterisk.  Thus, to assign the same sorption parameters to all components, only one line per model would be supplied, containing an asterisk in the first column.</p>
<p>Below is an example <cite>sorp</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">sorp</span> <span class="n">ltype</span>   <span class="n">a1l</span>     <span class="n">a2l</span>     <span class="n">bl</span>      <span class="n">vtype</span>   <span class="n">a1v</span>     <span class="n">a2v</span>     <span class="n">bv</span>
<span class="o">.</span><span class="n">model1</span>
<span class="n">CO3</span><span class="o">--</span>        <span class="n">freu</span>    <span class="mf">3.02</span>    <span class="mf">0.061</span>   <span class="mf">0.89</span>    <span class="o">*</span>       <span class="o">*</span>       <span class="o">*</span>       <span class="o">*</span>
<span class="n">C6H6</span> <span class="n">lin</span>     <span class="mf">0.89</span>    <span class="mf">3.16</span>    <span class="mf">0.2</span>     <span class="n">freu</span>    <span class="mf">0.35</span>    <span class="mf">4.519</span>   <span class="mf">0.688</span>
<span class="n">CO2</span>  <span class="n">mfreu</span>   <span class="mf">1.20</span>    <span class="mf">3.31</span>    <span class="mf">0.4</span>     <span class="n">con</span>     <span class="mf">0.01</span>    <span class="mf">2.01</span>    <span class="mf">0.61</span>
<span class="o">.</span><span class="n">model2</span>
<span class="n">CO3</span><span class="o">--</span>        <span class="n">con</span>     <span class="mf">0.001</span>   <span class="mf">2.02</span>    <span class="mf">0.88</span>    <span class="o">*</span>       <span class="o">*</span>       <span class="o">*</span>       <span class="o">*</span>
<span class="n">C6H6</span> <span class="n">con</span>     <span class="mf">1.20</span>    <span class="mf">0.58</span>    <span class="mf">1.12</span>    <span class="n">lin</span>     <span class="mf">2.2</span>     <span class="mf">2.043</span>   <span class="mf">2.7</span>
<span class="n">CO2</span>  <span class="n">mfreu</span>   <span class="mf">3.006</span>   <span class="mf">1.0</span>     <span class="mf">9.8</span>     <span class="n">mfreu</span>   <span class="mf">0.229</span>   <span class="mf">3.434</span>   <span class="mf">2.33</span>
</pre></div>
</div>
<p>In this example, C6H6 is modeled using the linear sorption isotherm model with a liquid α,,1,, parameter of 0.89 and a vapor α,,1,, parameter of 0.35 in adsorption model <cite>model1</cite>, and modeled using the conservative solute model with a liquid α,,1,, parameter of 1.20 and a vapor α,,1,, parameter of 2.2 in model <cite>model2</cite>.</p>
<p>If <code class="docutils literal notranslate"><span class="pre">sorp</span></code> is omitted, it is assumed to contain zeros for all values except the β parameters, which are assumed to be 1.</p>
</div>
<div class="section" id="cden">
<h2><code class="docutils literal notranslate"><span class="pre">cden</span></code><a class="headerlink" href="#cden" title="Permalink to this headline">¶</a></h2>
<p>The <cite>cden</cite> block allows the user to input the molecular weights of aqueous and aqueous Henry’s Law components and have the code adjust the density of the water according to the concentrations of these components.  It should not be used if <cite>trxn</cite> is preceded by a <cite>cden</cite> macro, as some values may be modified.  <cite>cden</cite> accepts several lines, each consisting of the name of a master species and its molecular weight, separated by spaces, tabs, or commas.</p>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">cden</span></code> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">cden</span>
<span class="n">HCO3</span><span class="o">-</span>        <span class="mi">61</span>
<span class="n">CO3</span><span class="o">--</span>        <span class="mi">60</span>
<span class="n">C6H6</span> <span class="mi">78</span>
</pre></div>
</div>
<p>In this example, HCO3- is defined to have a molecular weight of 61, CO3– to have a molecular mass of 60, and C6H6 to have a molecular mass of 78.</p>
<p>If database lookup is enabled (see the <cite>lookup</cite> block below), one of the lines in <cite>cden</cite> may consist of only an asterisk.  If such a line is provided, all components that <cite>lookup</cite> dynamically imports will be inserted into the <cite>cden</cite> block with the appropriate molecular weights.  These imported molecular weights can be overridden by explicitly listing the component and its desired molecular weight on a separate line.</p>
</div>
<div class="section" id="group">
<h2><code class="docutils literal notranslate"><span class="pre">group</span></code><a class="headerlink" href="#group" title="Permalink to this headline">¶</a></h2>
<p>The information in this block is only used if reactions are enabled (see <code class="docutils literal notranslate"><span class="pre">ctrl</span></code> above).  It is used to group aqueous components that take part in rapid kinetic reactions.  On the line below the keyword <cite>group</cite>, place one line for each group.  Each line should contain the names of all aqueous components present in that group, separated by spaces, tabs, or commas.</p>
<p>Below is an example of the <cite>group</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">group</span>
<span class="n">H</span>
<span class="n">C_a</span> <span class="n">Cl_a</span>
<span class="n">Na_a</span> <span class="n">Ca_a</span> <span class="n">Ca_a2</span>
<span class="n">U238</span> <span class="n">Th234</span>
</pre></div>
</div>
<p>In this example, H is in its own group, C_a and Cl_a are in a group together, Na_a, Ca_a, and Ca_a2 are grouped together, and U238 and Th234 are grouped together.</p>
</div>
<div class="section" id="cplx">
<h2><code class="docutils literal notranslate"><span class="pre">cplx</span></code><a class="headerlink" href="#cplx" title="Permalink to this headline">¶</a></h2>
<p>The <cite>cplx</cite> block allows the user to specify instantaneous reactions that form aqueous complexes from the aqueous master species and non-aqueous components specified in <cite>comp</cite>.  One equation is specified on each line, using a slightly modified version of the standard reaction format detailed in the <cite>rxn</cite> block below.  The left side of the reaction should contain only the name of the aqueous complex, without a number denoting its stoichiometry (the stoichiometry of the aqueous complex must be 1).  The right side should contain stoichiometry/compound pairs as specified by the standard format.  If a compound needs to be removed to make the aqueous complex, negate its stoichiometry.  After the equation, two comma-separated values in the format <cite>variable=value</cite> (not padded by spaces) are required:  <cite>ckeq</cite> for the equilibrium constant of that complex, and <cite>heq</cite> for the enthalpy of the formation of that complex.  On the same line as the keyword <cite>cplx</cite>, the keyword <cite>log10</cite> may be provided to denote that the values for all constants in the block will be given as the base-ten logarithm of the actual values.  Asterisks supplied in place of the equilibrium constant or enthalpy value always signify a zero, even if the <cite>log10</cite> keyword is specified.</p>
<p>If the keyword <cite>equi</cite> is supplied on the same line as the block name, the <cite>equi</cite> block below will be consulted to calculate equilibrium constants as functions of temperature.  In this case, the <cite>log10</cite> keyword and any values given to the right of the equations will be ignored.</p>
<p>If the database lookup option is enabled (see the <cite>lookup</cite> block below), complexes from the <cite>% CPLX</cite> section of the database file may be imported using a line that omits the right-hand side of the equation, and possibly the equilibrium information as well.  If such a line is encountered, the code will search for the named complex in the lookup database, and import the complex’s constituents, stoichiometry, and equilibrium-related constants.  The <cite>comp</cite> block will be automatically updated to include all components required for the imported complexes.  Note that the <cite>equi</cite> option cannot be used if database lookup is enabled, as this conflicts with the five-parameter fit used by PHREEQC databases for calculating the equilibrium constant as a function of temperature.  These equilibrium values can be overridden by placing explicit values for ckeq and heq on the same line.</p>
<p>Below is an example of <cite>cplx</cite>:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">cplx</span> <span class="n">log10</span>
<span class="n">CaHCO3</span><span class="o">+</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">Ca</span><span class="o">++</span> <span class="o">+</span> <span class="mi">1</span> <span class="n">HCO3</span><span class="o">-</span>                   <span class="n">ckeq</span><span class="o">=-</span><span class="mf">13.456</span><span class="p">,</span> <span class="n">heq</span><span class="o">=</span><span class="mi">0</span>
<span class="n">OH</span><span class="o">-</span> <span class="o">=</span> <span class="o">-</span><span class="mi">1</span> <span class="n">H</span><span class="o">+</span>                                  <span class="n">ckeq</span><span class="o">=-</span><span class="mi">14</span><span class="p">,</span> <span class="n">heq</span><span class="o">=</span><span class="mi">0</span>
<span class="n">Ca2UO2</span><span class="p">(</span><span class="n">CO3</span><span class="p">)</span><span class="mi">3</span> <span class="o">=</span> <span class="mi">2</span> <span class="n">Ca</span><span class="o">++</span> <span class="o">+</span> <span class="mi">2</span> <span class="n">UO2</span> <span class="o">+</span> <span class="mi">3</span> <span class="n">HCO3</span><span class="o">-</span> <span class="o">+</span> <span class="o">-</span><span class="mi">3</span> <span class="n">H</span><span class="o">+</span>      <span class="n">ckeq</span><span class="o">=-</span><span class="mf">20.26</span><span class="p">,</span> <span class="n">heq</span><span class="o">=*</span>
<span class="n">MgHPO4</span>
</pre></div>
</div>
<p>In this example, there are four complexes.  The complex Ca2UO2(CO3)3 is made by combining 2 Ca++, 2 UO2, and 3 HCO3- and removing 3 H+.  The information for the complex MgHPO4 is automatically imported from the database given in the <cite>lookup</cite> block.  All six equilibrium and enthalpy values are given as the base-ten logarithms of their actual values.  The equilibrium constant for CaHCO3+ is 1.0×10^-13.456^.  The enthalpies of CaHCO3+ and OH- are 1.0×10^0^ = 1, but the enthalpy of Ca2UO2(CO3)3 is zero.</p>
</div>
<div class="section" id="equi">
<h2><code class="docutils literal notranslate"><span class="pre">equi</span></code><a class="headerlink" href="#equi" title="Permalink to this headline">¶</a></h2>
<p>This block allows equilibrium constants for aqueous complexes from the <cite>cplx</cite> block above to vary with temperature.  If the keyword <cite>equi</cite> is provided in the <cite>cplx</cite> block above, <cite>trxn</cite> will not read equilibrium constants in <cite>cplx</cite>; instead, the <cite>equi</cite> block will be used to determine the constants as a function of temperature.</p>
<p>Every aqueous complex appearing in <cite>cplx</cite> must be included in <cite>equi</cite>.  For each aqueous complex, a line should be provided that contains only the name of the complex, prefixed by a single period.  On the lines that follow, the equilibrium constants should be specified by placing on each line a temperature in o C, a tab or space, and the equilibrium constant at that temperature.  The number of temperatures need not be the same for each complex.</p>
<p>Below is an example of the <cite>equi</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">equi</span>
<span class="o">.</span><span class="n">CaHCO3</span><span class="o">+</span>
<span class="mi">0</span>    <span class="mf">1e-10</span>
<span class="mi">10</span>   <span class="mf">2e-10</span>
<span class="mi">40</span>   <span class="mf">4e-9</span>
<span class="o">.</span><span class="n">OH</span><span class="o">-</span>
<span class="mi">0</span>    <span class="mf">1e-14</span>
<span class="o">.</span><span class="n">CaUO2</span><span class="p">(</span><span class="n">CO3</span><span class="p">)</span><span class="mi">3</span>
<span class="mi">20</span>   <span class="mf">1.3e-8</span>
<span class="mi">40</span>   <span class="mf">2.64e-6</span>
<span class="mi">60</span>   <span class="mf">8.7e-4</span>
<span class="mi">80</span>   <span class="mf">2.1e-1</span>
</pre></div>
</div>
<p>In this example, CaHCO3+ has an equilbirium constant of 1×10^-10^ at 0o C, 2×10^-10^ at 10o C, and 4×10^-9^ at 40o C.  OH- has an equilibrium constant of 1×10^-14^ at all temperatures.</p>
</div>
<div class="section" id="dist">
<h2><code class="docutils literal notranslate"><span class="pre">dist</span></code><a class="headerlink" href="#dist" title="Permalink to this headline">¶</a></h2>
<p>This block specifies distribution models that can be used for reaction types 1 and 2 to describe the distribution coefficient as a function of temperature.  For each distribution model, a line beginning with a period and then the name of the model (without a separating space) is given.  This is followed by a set of temperature/distribution coefficient pairs, one per line, for as many lines as desired.  FEHM will perform a piecewise linear interpolation between the values given to determine the value of the distribution coefficient for intervening temperatures.</p>
<p>Below is an example of the <cite>dist</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">dist</span>
<span class="o">.</span><span class="n">model1</span>
<span class="mi">0</span>    <span class="mi">1</span>
<span class="mi">10</span>   <span class="mi">2</span>
<span class="mi">20</span>   <span class="mi">4</span>
<span class="mi">30</span>   <span class="mi">8</span>
<span class="mi">40</span>   <span class="mi">16</span>
<span class="o">.</span><span class="n">model2</span>
<span class="mi">0</span>    <span class="mi">20</span>
<span class="mi">50</span>   <span class="mi">40</span>
<span class="mi">100</span>  <span class="mi">70</span>
</pre></div>
</div>
<p>There are two models in this example.  The first one, <cite>model1</cite>, has five data points at 10o C intervals from 0 to 40o C.  The distribution coefficient at 0o C is 1, the coefficient at 10o C is 2, and the coefficient at 20o C is 4.  Intermediate temperatures are linearly interpolated, so the distribution coefficient for <cite>model1</cite> at 25o C is 6.</p>
</div>
<div class="section" id="sol">
<h2><code class="docutils literal notranslate"><span class="pre">sol</span></code><a class="headerlink" href="#sol" title="Permalink to this headline">¶</a></h2>
<p>The <code class="docutils literal notranslate"><span class="pre">sol</span></code> block specifies solubility models that can be used for reaction types 7 and 8.  For each solubility model, a line beginning with a period, immediately followed by the model name, is given.  This is followed by temperature/solubility coefficient pairs, one per line, for as many lines as desired.  FEHM will perform a piecewise linear interpolation between the values given to determine the solubility coefficient for intervening temperatures.</p>
<p>Below is an example of the <cite>sol</cite> block:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">sol</span>
<span class="o">.</span><span class="n">model1</span>
<span class="mi">0</span>    <span class="mi">0</span>
<span class="mi">10</span>   <span class="mf">1e-16</span>
<span class="mi">20</span>   <span class="mf">1e-4</span>
<span class="mi">30</span>   <span class="mf">1e-2</span>
<span class="mi">50</span>   <span class="mi">1</span>
</pre></div>
</div>
<p>This example contains one model, <cite>model1</cite>.  In this model, the solubility coefficient changes from 0 to 1 over a range from 0o C to 50o C.  The solubility coefficient at 10o C is 1×10^-16^, and the coefficient at 20o C is 1×10^-4^.  Because FEHM linearly interpolates between successive values of the solubility coefficient, the coefficient at 40o C is 0.505.</p>
</div>
<div class="section" id="lookup">
<h2><code class="docutils literal notranslate"><span class="pre">lookup</span></code><a class="headerlink" href="#lookup" title="Permalink to this headline">¶</a></h2>
<p>This block enables the dynamic lookup process for mineral dissolution and aqueous complexation.  <cite>lookup</cite> uses a lookup database (generated by <cite>trxndb</cite> from a USGS PHREEQC geochemical database) to determine which reactions occur among a specific set of minerals and aqueous complexes.  The <cite>lookup</cite> block adds data to the <cite>comp</cite>, <cite>cden</cite>, <cite>group</cite>, <cite>cplx</cite>, and <cite>rxn</cite> blocks based on the information provided it.</p>
<p>In order to use <cite>lookup</cite>, a database in <cite>trxn</cite>’s standard format must be supplied.  This format is as follows:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span>Keyword `% MASTER`\\
Blank line(s)\\
Component 1 parameters:  component 1 name, component 1 master species name, component 1 molar mass\\
Component 2 parameters...\\
Blank line(s)\\
Keyword `% CPLX`\\
Blank line(s)\\
Complex 1 name\\
Complex 1 products:  product 1 stoichiometry, product 1 name; product 2 stoichiometry, product 2 name...\\
Equilibrium constant for this complex\\
Enthalpy of reaction for this complex\\
Up to five values defining temperature dependence of equilibrium constant:  A,,1,,, A,,2,,, A,,3,,, A,,4,,, A,,5,,, where log,,10,, K = A,,1,, + A,,2,, · T + A,,3,, / T + A,,4,, · log,,10,, T + A,,5,, / T^2^.\\
Blank line(s)\\
Complex 2 parameters...\\
Blank line(s)\\
Keyword `% MIN`\\
Blank line(s)\\
Mineral 1 name, mineral 1 formula\\
Mineral 1 products\\
Equilibrium constant for this mineral\\
Enthalpy of reaction for this mineral\\
Up to five values defining temperature dependence of equilibrium constant\\
Blank line(s)\\
Mineral 2 parameters...\\
Blank line(s)\\
Keyword `% END`\\
</pre></div>
</div>
<p>Comments can be included anywhere in the input file by using a pound sign (“#”).</p>
<p>Below is a brief example of a database in this form.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="o">%</span> <span class="n">MASTER</span>

<span class="n">E</span>    <span class="n">e</span><span class="o">-</span>      <span class="mi">0</span>
<span class="n">H</span>    <span class="n">H</span><span class="o">+</span>      <span class="mf">1.007942</span>
<span class="n">Mn</span><span class="p">(</span><span class="o">+</span><span class="mi">2</span><span class="p">)</span>       <span class="n">Mn</span><span class="o">++</span>    <span class="mf">54.938045</span>
<span class="n">Mn</span><span class="p">(</span><span class="o">+</span><span class="mi">3</span><span class="p">)</span>       <span class="n">Mn</span><span class="o">+++</span>   <span class="mf">54.938045</span>
<span class="n">F</span>    <span class="n">F</span><span class="o">-</span>      <span class="mf">18.9984032</span>
<span class="n">Al</span>   <span class="n">AlOH</span><span class="o">++</span>  <span class="mf">26.9815386</span>
<span class="n">Si</span>   <span class="n">H4SiO4</span>  <span class="mf">28.08554</span>
<span class="n">Mg</span>   <span class="n">Mg</span><span class="o">++</span>    <span class="mf">24.305</span>
<span class="n">Ca</span>   <span class="n">Ca</span><span class="o">++</span>    <span class="mf">40.0782</span>
<span class="n">S</span>    <span class="n">HS</span><span class="o">-</span>     <span class="mf">32.0652</span>
<span class="n">Fe</span>   <span class="n">Fe</span><span class="o">++</span>    <span class="mf">55.845</span>

<span class="o">%</span> <span class="n">CPLX</span>

<span class="n">Al</span><span class="o">+++</span>
<span class="mi">1</span> <span class="n">AlOH</span><span class="o">++</span> <span class="mi">1</span> <span class="n">H</span><span class="o">+</span>
<span class="o">-</span><span class="mf">5.00</span>
<span class="mf">11.49</span>
 <span class="o">-</span><span class="mf">38.253</span> <span class="mf">0.0</span> <span class="o">-</span><span class="mf">656.27</span> <span class="mf">14.327</span>

<span class="n">Mn</span><span class="o">++</span>
<span class="mi">1</span> <span class="n">e</span><span class="o">-</span> <span class="mi">1</span> <span class="n">Mn</span><span class="o">+++</span>
<span class="o">-</span><span class="mf">25.510</span>
<span class="mf">25.800</span>


<span class="n">MnF</span><span class="o">+</span>
<span class="mi">1</span> <span class="n">Mn</span><span class="o">++</span> <span class="mi">1</span> <span class="n">F</span><span class="o">-</span>
<span class="mf">0.840</span>

<span class="o">%</span> <span class="n">IMM</span>

<span class="n">Sepiolite</span> <span class="n">Mg2Si3O7</span><span class="o">.</span><span class="mi">5</span><span class="n">OH</span><span class="p">:</span><span class="mi">3</span><span class="n">H2O</span>
<span class="mi">3</span> <span class="n">H4SiO4</span> <span class="o">-</span><span class="mi">4</span> <span class="n">H</span><span class="o">+</span> <span class="mi">2</span> <span class="n">Mg</span><span class="o">++</span>
<span class="mf">15.760</span>
<span class="o">-</span><span class="mf">10.700</span>

<span class="n">Fluorite</span> <span class="n">CaF2</span>
<span class="mi">1</span> <span class="n">Ca</span><span class="o">++</span> <span class="mi">2</span> <span class="n">F</span><span class="o">-</span>
<span class="o">-</span><span class="mf">10.600</span>
<span class="mf">4.690</span>
<span class="mf">66.348</span> <span class="mf">0.0</span> <span class="o">-</span><span class="mf">4298.2</span> <span class="o">-</span><span class="mf">25.271</span>


<span class="n">Pyrite</span> <span class="n">FeS2</span>
<span class="mi">2</span> <span class="n">HS</span><span class="o">-</span> <span class="o">-</span><span class="mi">2</span> <span class="n">e</span><span class="o">-</span> <span class="mi">1</span> <span class="n">Fe</span><span class="o">++</span> <span class="o">-</span><span class="mi">2</span> <span class="n">H</span><span class="o">+</span>
<span class="o">-</span><span class="mf">18.479</span>
<span class="mf">11.300</span>

<span class="o">%</span> <span class="n">END</span>
</pre></div>
</div>
<p>In this example, there are three complexes and three minerals.  The first complex is Al+++, which is formed by the master species AlOH++ and H+.  Its equilibrium constant is -5.00, its enthalpy of formation is 11.49, and four of five possible temperature-dependence parameters are supplied.  For the complex MnF+, no enthalpy of formation or temperature-dependence parameters are supplied.</p>
<p>A conversion script, <code class="docutils literal notranslate"><span class="pre">trxndb</span></code>, is available to automate conversion of certain USGS PHREEQC input files to the appropriate format.  The converter script understands a limited subset of the complete syntax used in PHREEQC input files; if it gives improper results or errors, ensure that the input file is in a consistent format and that the keywords <cite>log_k</cite>, <cite>delta_h</cite>, and <cite>-analytic</cite> are used rather than their shortened alternatives.  An example file that can be converted flawlessly by <cite>trxndb</cite> can be found at <cite>/scratch/nts/ms/trxn/geochem/phreeqc.dat</cite>.  The converter script is written in Perl and located at <cite>/scratch/nts/ms/trxn/geochem/trxndb</cite>.  The script should be called with the name of the PHREEQC input file as the only argument.  It will create a file in the same directory with the same name and extension <cite>trxd</cite> containing the <cite>trxn</cite>-compatible database.  For example:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span>$ ls
phreeqc.dat
$ /scratch/nts/ms/trxn/geochem/trxndb phreeqc.dat
45 master species, 187 solution equations read
58 mineral equations read
45 master species, 180 solution equations written
57 mineral equations written
$ ls
phreeqc.dat  phreeqc.trxd
$
</pre></div>
</div>
<p>The <code class="docutils literal notranslate"><span class="pre">lookup</span></code> block must be the first block in the <code class="docutils literal notranslate"><span class="pre">trxn</span></code> macro (excepting the <code class="docutils literal notranslate"><span class="pre">ctrl</span></code> block if it is used) and consists of one line, which contains the block name <cite>lookup</cite> and the full path to the database file.  The <cite>lookup</cite> block should be followed by a blank or commented line.</p>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">lookup</span></code> block.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">lookup</span> <span class="o">/</span><span class="n">scratch</span><span class="o">/</span><span class="n">nts</span><span class="o">/</span><span class="n">ms</span><span class="o">/</span><span class="n">trxn</span><span class="o">/</span><span class="n">geochem</span><span class="o">/</span><span class="n">phreeqc</span><span class="o">.</span><span class="n">trxd</span>
</pre></div>
</div>
<p>In this example, the database is located at <code class="docutils literal notranslate"><span class="pre">/scratch/nts/ms/trxn/geochem/phreeqc.trxd</span></code>.</p>
<p>Once <code class="docutils literal notranslate"><span class="pre">lookup</span></code> has been enabled, the <code class="docutils literal notranslate"><span class="pre">cden</span></code>, <code class="docutils literal notranslate"><span class="pre">cplx</span></code>, and <code class="docutils literal notranslate"><span class="pre">rxn</span></code> blocks may be modified by the user to utilize the information from the lookup database.  Please see the sections for these blocks for the appropriate syntax to take advantage of this information.  If the <cite>debug</cite> option is provided in the <cite>ctrl</cite> block, the final version of each block will be printed to the output (<cite>.out</cite>) file, which may aid in debugging if the results of the run are not as expected.</p>
</div>
<div class="section" id="rxn">
<h2><code class="docutils literal notranslate"><span class="pre">rxn</span></code><a class="headerlink" href="#rxn" title="Permalink to this headline">¶</a></h2>
<p>The <code class="docutils literal notranslate"><span class="pre">rxn</span></code> blocks are used to model kinetic reactions between simulated compounds.  Seven types of reactions my be used, each with its own input parameters and input format.  <cite>rxn</cite> is intended to be specified multiple times, once for each reaction that is taking place.  For each reaction, the first line of the block contains the keyword <cite>rxn</cite> and the number representing the type of the reaction (see below).  The next several lines are used for the parameters unique to that reaction type, which are detailed below.  The end of the <cite>rxn</cite> block is signaled by a blank line.</p>
<p>Most of these reactions take one line of input in the standard reaction format, which is used to specify reactants, products, and stoichiometries simultaneously.  In this format, each reactant or product is specified by a number denoting the stoichiometry of that compound, a space, and then the name of the compound as given in <cite>comp</cite>.  Compounds are separated from each other by a plus sign padded on either side with spaces (” + “).  The products and reactants are separated from each other by a token containing an equals sign (“=”).  Optionally, directional indicators may be added to the equals sign to indicate the direction of the reaction (e.g., “=&gt;”).  The reactants must be placed on the left side of the equals sign, and the products on the right side.  For example:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="mi">6</span> <span class="n">HCl</span> <span class="o">+</span> <span class="mi">2</span> <span class="n">Al_s</span> <span class="o">=&gt;</span> <span class="mi">2</span> <span class="n">AlCl3</span> <span class="o">+</span> <span class="mi">3</span> <span class="n">H2</span>
</pre></div>
</div>
<p>The reaction types and their parameters are as follows.  A depthier description of the mechanics of each type of reaction can be found at the end of the <cite>rxn</cite> section of the FEHM User’s Manual.</p>
<table border="1" class="docutils">
<colgroup>
<col width="1%" />
<col width="2%" />
<col width="97%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Number</th>
<th class="head">Type</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>Linear kinetic reaction</td>
<td>This reaction accepts one aqueous reactant and one solid product with 1:1 stoichiometry.  The first parameter line for this reaction is a reaction in standard format, without the stoichiometric coefficients.  The second parameter line in a list of comma-separated name/value pairs, where the acceptable names are <cite>rate</cite> to specify the rate of the reaction and <cite>distcoef</cite> to specify the distribution coefficient.  The value for <cite>distcoef</cite> my be a real number of the coefficient is constant, or the name of a model specified in <cite>dist</cite> (without the beginning period).</td>
</tr>
<tr class="row-odd"><td>2</td>
<td>Langmuir kinetic reaction</td>
<td>This reaction’s parameter input is identical to the input for reaction type 1 except for an extra available parameter for the name/value pairs.  This parameter is <cite>maxconc</cite>, used to set the maximum concentration that can sorb.</td>
</tr>
<tr class="row-even"><td>3</td>
<td>General kinetic reaction</td>
<td>This reaction accepts its first line of of input as a generic reaction in the standard format described above.  This is followed by a line of name/value pairs, for which the name can be <cite>forward</cite> to set the forward rate constant for the reaction or <cite>reverse</cite> to set the reverse rate constant.</td>
</tr>
<tr class="row-odd"><td>4</td>
<td>Dual Monod biodegradation reaction</td>
<td>The first three to six parameter lines of this reaction consist of a name, a colon, and one or more component/complex names.  The name <cite>substrate</cite> accepts a single immobile component name that is to be the substrate that is degraded.  The name <cite>electronacceptor</cite> accepts the stoichiometry and name (separated by a space) of the aqueous complex that is the electron acceptor for the reaction.  The name <cite>biomass</cite> accepts the stoichiometry and name of the solid component that is the biomass produced by the reaction.  The name <cite>reactants</cite> is optional and accepts the stoichiometries and names of any extra reactants that are participating in the reaction.  Likewise, <cite>products</cite> is optional and accepts the names of any extra products of the reaction.  Only a total of three additional reactants and products can be specified.  If only certain forms of the substrate are biodegradable, those can be listed with the <cite>biodegradable</cite> name.\The last parameter line contains a list of name/value pairs as follows:  <cite>ks</cite> - the half-maximum-rate concentration of the substrate; <cite>ka</cite> - the half-maximum-rate concentration of the electron acceptor; <cite>decay</cite> - the microbial decay rate (1/hr); <cite>phthreshold</cite> - the pH threshold above which biodegradation will not occur; <cite>qm</cite> - the maximum rate of substrate utilization (mol/kg biomass/hr); <cite>yield</cite> - the microbial yield coefficient (kg biomass/mol substrate); <cite>xminit</cite> - the minimum concentration of biomass (mol/kg rock).</td>
</tr>
<tr class="row-even"><td>5</td>
<td>Radioactive decay reaction</td>
<td>This reaction accepts one aqueous component reactant and one aqueous component product with 1:1 stoichiometry.  This first parameter line is a reaction in standard format without stoichiometric coefficients.  The reactant is the component that is decaying, and the product is the decay product.  If the decay product is not being modeled in the simulation, an asterisk may be given in place of the product name.  The second parameter line contains a single name/value pair.  The name is <cite>halflife</cite>, and the value is the half-life of the reaction in years.</td>
</tr>
<tr class="row-odd"><td>7 and 8</td>
<td>Precipitation/dissolution reaction</td>
<td>For reaction type 7, the rates are based on the total aqueous concentration of the components; whereas for reaction type 8, the rates are based on the free-ion concentrations alone.  This reaction accepts one solid component and any number of aqueous master species.  The first line of input for this reaction should be a reaction in standard format, containing only a solid component on one side, and at least one master species on the other.  The second line should contain the following name/value pairs:  <cite>solubility</cite> - the solubility product (either a real number or the name of a solubility model specified in <cite>sol</cite>); <cite>rate</cite> - the reaction rate (mol/m^2^/s); <cite>sarea</cite> - the surface area of the mineral (m^2^/m^3^ rock); <cite>molecularweight</cite> - the molecular weight of the mineral; <cite>density</cite> - the density of the mineral.  <cite>molecularweight</cite> and <cite>density</cite> should only be provided for reaction type 8.\For reaction types 7 and 8 only, if the database lookup option is enabled (see the <cite>lookup</cite> block above), an alternate form of the reaction can be input.  On the first line of the reaction, provide only the name of a mineral defined in the <cite>% IMM</cite> section of the lookup database.  The information on the reactants and products, as well as the solubility constant, are imported from the database, and any components that are necessary for the reaction are dynamically added to the <cite>comp</cite> block.  However, the <cite>rate</cite> and <cite>sarea</cite> parameters (and <cite>molecularweight</cite> and <cite>density</cite> parameters for reaction type 8) are not contained in the database, and must still be set by the user on the second line of the reaction.</td>
</tr>
</tbody>
</table>
<p>Please note that in all name/value pairs, the name and value must be separated by only an equals sign that is not padded by spaces.</p>
<p>Here is an example of each reaction type.  Two versions of reaction type 5 are given to demonstrate the optional unsimulated daughter species, and two versions of reaction type 7 are provided to demonstrate the database lookup option.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">rxn</span> <span class="mi">1</span>
<span class="n">Ca_a</span> <span class="o">&lt;=&gt;</span> <span class="n">Ca_s</span>
<span class="n">rate</span><span class="o">=</span><span class="mi">10</span><span class="p">,</span> <span class="n">distcoef</span><span class="o">=</span><span class="mi">2</span>

<span class="n">rxn</span> <span class="mi">2</span>
<span class="n">Ca_a</span> <span class="o">&lt;=&gt;</span> <span class="n">Ca_s</span>
<span class="n">rate</span><span class="o">=</span><span class="mi">10</span><span class="p">,</span> <span class="n">maxconc</span><span class="o">=</span><span class="mi">2</span><span class="p">,</span> <span class="n">distcoef</span><span class="o">=</span><span class="n">model2</span>

<span class="n">rxn</span> <span class="mi">3</span>
<span class="mi">3</span> <span class="n">H</span> <span class="o">+</span> <span class="mi">2</span> <span class="n">Ca</span> <span class="o">+</span> <span class="mi">5</span> <span class="n">U238</span> <span class="o">&lt;=&gt;</span> <span class="n">Cl</span> <span class="o">+</span> <span class="mi">2</span> <span class="n">Cl2</span> <span class="o">+</span> <span class="mi">3</span> <span class="n">C_a</span>
<span class="n">forward</span><span class="o">=</span><span class="mi">2</span><span class="p">,</span> <span class="n">reverse</span><span class="o">=</span><span class="mf">1.5</span>

<span class="n">rxn</span> <span class="mi">4</span>
<span class="n">substrate</span><span class="p">:</span>  <span class="n">U238</span>
<span class="n">electronacceptor</span><span class="p">:</span>  <span class="mi">3</span> <span class="n">H</span>
<span class="n">biomass</span><span class="p">:</span>  <span class="n">Ca</span>
<span class="n">reactants</span><span class="p">:</span>  <span class="mi">2</span> <span class="n">Cl</span>
<span class="n">products</span><span class="p">:</span>  <span class="mi">1</span> <span class="n">Na</span><span class="p">,</span> <span class="mi">5</span> <span class="n">Th234</span>
<span class="n">biodegradable</span><span class="p">:</span>  <span class="n">UO2</span>
<span class="n">ks</span><span class="o">=</span><span class="mf">1.2</span><span class="p">,</span> <span class="n">ka</span><span class="o">=</span><span class="mf">1.35</span><span class="p">,</span> <span class="n">decay</span><span class="o">=</span><span class="mf">0.69</span><span class="p">,</span> <span class="n">phthreshold</span><span class="o">=</span><span class="mf">8.2</span><span class="p">,</span> <span class="n">qm</span><span class="o">=</span><span class="mf">0.20</span><span class="p">,</span> <span class="k">yield</span><span class="o">=</span><span class="mf">1.2</span><span class="p">,</span> <span class="n">xminit</span><span class="o">=</span><span class="mf">0.067</span>

<span class="n">rxn</span> <span class="mi">5</span>
<span class="n">U238</span> <span class="o">=&gt;</span> <span class="n">Th234</span>
<span class="n">halflife</span><span class="o">=</span><span class="mi">20</span>

<span class="n">rxn</span> <span class="mi">5</span>
<span class="n">U234</span> <span class="o">=&gt;</span> <span class="o">*</span>
<span class="n">halflife</span><span class="o">=</span><span class="mi">10</span>

<span class="n">rxn</span> <span class="mi">7</span>
<span class="n">NaCl</span> <span class="o">&lt;=&gt;</span> <span class="n">Na</span> <span class="o">+</span> <span class="n">Cl</span>
<span class="n">solubility</span><span class="o">=</span><span class="mf">0.0231</span><span class="p">,</span> <span class="n">rate</span><span class="o">=</span><span class="mf">1.02</span><span class="p">,</span> <span class="n">sarea</span><span class="o">=</span><span class="mf">2.2</span>

<span class="n">rxn</span> <span class="mi">7</span>
<span class="n">Quartz</span>
<span class="n">rate</span><span class="o">=</span><span class="mf">0.2</span><span class="p">,</span> <span class="n">sarea</span><span class="o">=</span><span class="mi">1</span>

<span class="n">rxn</span> <span class="mi">8</span>
<span class="n">CaCl2</span> <span class="o">&lt;=&gt;</span> <span class="n">Ca</span> <span class="o">+</span> <span class="mi">2</span> <span class="n">Cl</span>
<span class="n">solubility</span><span class="o">=</span><span class="n">model1</span><span class="p">,</span> <span class="n">rate</span><span class="o">=</span><span class="mf">0.2</span><span class="p">,</span> <span class="n">sarea</span><span class="o">=</span><span class="mi">5</span><span class="p">,</span> <span class="n">molecularweight</span><span class="o">=</span><span class="mi">60</span><span class="p">,</span> <span class="n">density</span><span class="o">=</span><span class="mf">5.25</span>
</pre></div>
</div>
</div>
<div class="section" id="assign">
<h2><code class="docutils literal notranslate"><span class="pre">assign</span></code><a class="headerlink" href="#assign" title="Permalink to this headline">¶</a></h2>
<p>This block allows the user to assign the parameters stored in models in the above blocks to the zones defined in the <cite>zone</cite> macro.  The <cite>assign</cite> block also allows assignment of some other parameters that are specific to zones.</p>
<p>Zone numbers run down the side of the <cite>assign</cite> block, and parameters run across the top.  The following parameters may be supplied; all are optional:</p>
<table border="1" class="docutils">
<colgroup>
<col width="1%" />
<col width="87%" />
<col width="12%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Option</th>
<th class="head">Description</th>
<th class="head">Default Value</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><cite>water</cite></td>
<td>Sets the initial water type filling the zone.  This parameter must be the name of a water type defined in the <cite>water</cite> block.</td>
<td>Pure water</td>
</tr>
<tr class="row-odd"><td><cite>boun</cite></td>
<td>Sets the concentrations of species coming from inflow nodes in this zone.  This parameter may consist of a water type defined in <cite>water</cite> and/or a gas type defined in <cite>gas</cite>.  If both a water type and a gas type are flowing in, they should be separated by a period (not padded by spaces).  If there is no inflow in this zone, give an asterisk for this parameter.  An asterisk can also be used to specify inflow of pure water.</td>
<td>Pure water</td>
</tr>
<tr class="row-even"><td><cite>time</cite></td>
<td>Sets the time range during which the inflow nodes in this zone are injecting, in days.  A lone zero (“0”) gives no inflow.  An asterisk provides no inflow if the entry in the <cite>boun</cite> column is “*” or if the <cite>boun</cite> column is missing, but provides inflow over the entire simulation otherwise.  A single number other than zero will give inflow for one day starting at the specified day.  A range separated by a greater-than sign (“&gt;”) not padded by spaces will run injection from the first number specified to the second number.  If the number before or after the greater-than is ommitted, it will default to the beginning or end of the simulation, respectively.  Thus, “&gt;30” will run injection from time 0 through time 30, “30&gt;” will run injection from time 30 to the end of the simulation, and “&gt;” will run injection for the entire simulation.</td>
<td>Inflow for the entire simulation if valid water/gas types are specified in the <cite>boun</cite> column, no inflow otherwise</td>
</tr>
<tr class="row-odd"><td><cite>rock</cite></td>
<td>Sets the composition of the rock in this zone.  This parameter must be the name of a rock type defined in <cite>rock</cite>.</td>
<td>Rock contains no relevant species</td>
</tr>
<tr class="row-even"><td><cite>gas</cite></td>
<td>Sets the initial composition of the gas in this zone.  This parameter must be the name of a gas type defined in <cite>gas</cite>.  If there is no gas in this zone, give an asterisk for this parameter.</td>
<td>No gas</td>
</tr>
<tr class="row-odd"><td><cite>disp</cite></td>
<td>Sets the dispersivity constants for this zone.  This parameter must be the name of a dispersivity model defined in <cite>disp</cite>.</td>
<td>No dispersivity</td>
</tr>
<tr class="row-even"><td><cite>sorp</cite></td>
<td>Sets adsorption parameters for this zone.  This parameter must be the name of an adsorption model defined in <cite>sorp</cite>.  While the models in <cite>sorp</cite> are defined by starting the line with a period, do not include the period in this parameter.</td>
<td>No adsorption for any components</td>
</tr>
<tr class="row-odd"><td><cite>tpor</cite></td>
<td>Sets the optional tracer porosity for this zone.  This parameter must be a real number from 0 to 1.</td>
<td>0.32</td>
</tr>
<tr class="row-even"><td><cite>opt</cite></td>
<td>Sets miscellaneous options per-zone.  Options should be separated by periods with no spaces.  The available options are:  <cite>const</cite> causes the concentrations of solutes in the inflow water to be held constant at nodes in the current zone; <cite>accum</cite> enables solute accumulation in the current zone.  Note that <cite>const</cite> and <cite>accum</cite> are mutually exclusive.  An asterisk can be used to specify no options.</td>
<td>No options</td>
</tr>
</tbody>
</table>
<p>An asterisk may be given for any of the parameters above for a given zone, in which case the default values given above are used.  If an entire column is omitted, that parameter will be given the default values shown above for every zone.  If a zone is omitted, it will receive default values for every column.  Nodes that are not in any zone at the time when <cite>trxn</cite> is called will also receive default values.</p>
<p>Below is an example of the <code class="docutils literal notranslate"><span class="pre">assign</span></code> block.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">assign</span>       <span class="n">water</span>   <span class="n">rock</span>    <span class="n">gas</span>     <span class="n">boun</span>    <span class="n">time</span>    <span class="n">disp</span>    <span class="n">sorp</span>    <span class="n">tpor</span>
<span class="mi">1</span>    <span class="n">wt3</span>     <span class="n">granite</span> <span class="o">*</span>       <span class="n">vt1</span><span class="o">.</span><span class="n">wt1</span> <span class="mi">20</span><span class="o">&gt;</span><span class="mi">30</span>   <span class="n">model1</span>  <span class="n">model1</span>  <span class="mf">0.28</span>
<span class="mi">2</span>    <span class="n">wt2</span>     <span class="n">clay</span>    <span class="n">vt2</span>     <span class="o">*</span>       <span class="o">*</span>       <span class="n">model2</span>  <span class="n">model2</span>  <span class="mf">0.69</span>
<span class="mi">3</span>    <span class="n">wt2</span>     <span class="n">granite</span> <span class="n">vt2</span>     <span class="o">*</span>       <span class="o">*</span>       <span class="n">model1</span>  <span class="n">model2</span>  <span class="mf">0.32</span>
</pre></div>
</div>
<p>In this example, there are three zones.  The first zone is initially filled with water of type <cite>wt3</cite> and no gases, and water of type <cite>wt1</cite> and gas of type <cite>vt1</cite> are flowing into it starting at time 20 days and ending at time 30 days.  Its rock is of type <cite>granite</cite>, it uses the dispersivity parameters defined in <cite>model1</cite> and the adsorption parameters defined in <cite>model1</cite>, and has a tracer porosity of 0.28.  Note also that zones 2 and 3 have no inflow.</p>
<p>After specifying all applicable blocks, use the <cite>end trxn</cite> keyword to end the reading of the <cite>trxn</cite> macro.</p>
</div>
<div class="section" id="example">
<h2>Example<a class="headerlink" href="#example" title="Permalink to this headline">¶</a></h2>
<p>Below is a complete, commented working example <cite>trxn</cite> macro, with its accompnaying zone macro for reference.   This example was taken from the <cite>multi_solute</cite> test case in FEHM’s standard verification suite.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">zone</span>
<span class="n">default</span>
<span class="nb">all</span>
<span class="n">bound</span>
<span class="n">nnum</span>
<span class="mi">2</span>    <span class="mi">1</span> <span class="mi">102</span>

<span class="n">trxn</span>

<span class="c1"># In the zone macro above, all nodes are placed into the zone named &quot;default&quot;,</span>
<span class="c1"># except for nodes 1 and 202, which are placed in the zone named &quot;bound&quot;.</span>
<span class="c1"># Zone &quot;bound&quot; will be used for inflow.</span>

<span class="n">ctrl</span> <span class="n">rxnon</span> <span class="c1"># Enable reactions.</span>

<span class="c1"># Include header information from trac.</span>
<span class="n">header</span>
<span class="mf">1.</span><span class="n">d</span><span class="o">-</span><span class="mi">80</span>       <span class="mf">1.0</span>     <span class="mf">1.e-6</span>   <span class="mf">0.5</span>
<span class="mf">1.</span>   <span class="mi">2000</span>    <span class="mf">1.0</span>     <span class="mi">2000</span>
<span class="mi">5</span>    <span class="mf">5.0</span>     <span class="mf">1.e-6</span>   <span class="mf">2.8333e-3</span>       <span class="mi">0</span>
<span class="n">iskip</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">rsdmax</span><span class="o">=</span><span class="mf">1e-10</span>

<span class="c1"># There are six components in this simulation:  aqueous cobalt, iron, and EDTA,</span>
<span class="c1"># and solid Co-EDTA, Fe-EDTA, and cobalt.</span>
<span class="n">comp</span>         <span class="n">master</span>
<span class="n">a</span> <span class="n">Cobalt</span>     <span class="n">Cobalt_a</span>
<span class="n">a</span> <span class="n">Iron</span>               <span class="n">Iron_a</span>
<span class="n">a</span> <span class="n">EDTA</span>               <span class="n">EDTA_a</span>
<span class="n">s</span> <span class="n">Co</span><span class="o">-</span><span class="n">EDTA_s</span>  <span class="o">*</span>
<span class="n">s</span> <span class="n">Fe</span><span class="o">-</span><span class="n">EDTA_s</span>  <span class="o">*</span>
<span class="n">s</span> <span class="n">Cobalt_s</span>   <span class="o">*</span>

<span class="c1"># There is only one type of water, called &quot;inflow&quot;, which contains 3.16e-5 M</span>
<span class="c1"># aqueous cobalt, 1e-13 M aqueous iron, and 3.16e-5 M aqueous EDTA.</span>
<span class="n">water</span>        <span class="n">Cobalt</span>  <span class="n">Iron</span>    <span class="n">EDTA</span>
<span class="n">inflow</span>       <span class="mf">3.16e-5</span> <span class="mf">1e-13</span>   <span class="mf">3.16e-5</span>

<span class="c1"># There is one sorption model.  It models all components with a linear sorption</span>
<span class="c1"># isotherm, using alpha-1 and alpha-2 parameters of zero, and a beta parameter</span>
<span class="c1"># of 1.</span>
<span class="n">sorp</span> <span class="n">ltype</span>   <span class="n">a1l</span>     <span class="n">a2l</span>     <span class="n">bl</span>
<span class="o">.</span><span class="n">smod</span>
<span class="n">Cobalt</span>       <span class="n">lin</span>     <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>
<span class="n">Iron</span> <span class="n">lin</span>     <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>
<span class="n">EDTA</span> <span class="n">lin</span>     <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>

<span class="c1"># We assign the liquid diffusion coefficient for the simulation to be 1e-9.</span>
<span class="n">diff</span> <span class="n">l</span><span class="o">=</span><span class="mf">1e-9</span>

<span class="c1"># There is one model for dispersivity, &quot;dmod&quot;.  It sets the dispersivity to</span>
<span class="c1"># 0.05 in the X direction, and 1e-34 in the Y and Z directions.</span>
<span class="n">disp</span> <span class="n">lx</span>      <span class="n">ly</span>      <span class="n">lz</span>
<span class="n">dmod</span> <span class="mf">0.05</span>    <span class="mf">1e-34</span>   <span class="mf">1e-34</span>

<span class="c1"># There are two groups for the coupled solver.  One contains cobalt and EDTA,</span>
<span class="c1"># and the other contains iron.</span>
<span class="n">group</span>
<span class="n">Cobalt</span> <span class="n">EDTA</span>
<span class="n">Iron</span>

<span class="c1"># Here, two aqueous complexes are defined:  Co-EDTA and Fe-EDTA.  Co-EDTA is</span>
<span class="c1"># composed of one aqueous cobalt and one aqueous EDTA; Fe-EDTA is composed of</span>
<span class="c1"># one aqueous iron and one aqueous EDTA.  Both complexes have enthalpy zero;</span>
<span class="c1"># Co-EDTA has equilibrium constant 1e18 and Fe-EDTA has equilibrium constant</span>
<span class="c1"># 6.31e27.</span>
<span class="n">cplx</span>
<span class="n">Co</span><span class="o">-</span><span class="n">EDTA_a</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">Cobalt_a</span> <span class="o">+</span> <span class="mi">1</span> <span class="n">EDTA_a</span>    <span class="n">ckeq</span><span class="o">=</span><span class="mf">1e18</span><span class="p">,</span> <span class="n">heq</span><span class="o">=</span><span class="mi">0</span>
<span class="n">Fe</span><span class="o">-</span><span class="n">EDTA_a</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">Iron_a</span> <span class="o">+</span> <span class="mi">1</span> <span class="n">EDTA_a</span>              <span class="n">ckeq</span><span class="o">=</span><span class="mf">6.31e27</span><span class="p">,</span> <span class="n">heq</span><span class="o">=</span><span class="mi">0</span>

<span class="c1"># There are four reactions taking place in this simulation.</span>

<span class="c1"># Reaction 1 is a linear kinetic reaction describing the dissolution of</span>
<span class="c1"># Co-EDTA.  The distribution coefficient is 0.533, and the rate of reaction is</span>
<span class="c1"># 1.</span>
<span class="n">rxn</span> <span class="mi">1</span>
<span class="n">Co</span><span class="o">-</span><span class="n">EDTA_a</span> <span class="o">&lt;=&gt;</span> <span class="n">Co</span><span class="o">-</span><span class="n">EDTA_s</span>
<span class="n">distcoef</span><span class="o">=</span><span class="mf">0.533</span><span class="p">,</span> <span class="n">rate</span><span class="o">=</span><span class="mi">1</span>

<span class="c1"># Reaction 2, also linear kinetic, describes the dissolution of cobalt.  The</span>
<span class="c1"># distribution coefficient is 5.07, and the rate of reaction is again 1.</span>
<span class="n">rxn</span> <span class="mi">1</span>
<span class="n">Cobalt_a</span> <span class="o">&lt;=&gt;</span> <span class="n">Cobalt_s</span>
<span class="n">distcoef</span><span class="o">=</span><span class="mf">5.07</span><span class="p">,</span> <span class="n">rate</span><span class="o">=</span><span class="mi">1</span>

<span class="c1"># Reaction 3 is also linear kinetic and describes the dissolution of Fe-EDTA.</span>
<span class="c1"># The distribution coefficient here is 0.427, and the rate of reaction is 1.</span>
<span class="n">rxn</span> <span class="mi">1</span>
<span class="n">Fe</span><span class="o">-</span><span class="n">EDTA_a</span> <span class="o">&lt;=&gt;</span> <span class="n">Fe</span><span class="o">-</span><span class="n">EDTA_s</span>
<span class="n">distcoef</span><span class="o">=</span><span class="mf">0.427</span><span class="p">,</span> <span class="n">rate</span><span class="o">=</span><span class="mi">1</span>

<span class="c1"># Reaction 4 is a general kinetic reaction that describes the complexation of</span>
<span class="c1"># solid cobalt and Fe-EDTA to form Co-EDTA.  The forward rate constant is</span>
<span class="c1"># 1.26e-2, and the reaction never occurs in reverse (the reverse rate constant</span>
<span class="c1"># is zero).</span>
<span class="n">rxn</span> <span class="mi">3</span>
<span class="n">Co</span><span class="o">-</span><span class="n">EDTA_s</span> <span class="o">=</span> <span class="n">Fe</span><span class="o">-</span><span class="n">EDTA_s</span> <span class="o">+</span> <span class="n">Cobalt_s</span>
<span class="k">for</span><span class="o">=</span><span class="mf">1.26e-2</span><span class="p">,</span> <span class="n">rev</span><span class="o">=</span><span class="mi">0</span>

<span class="c1"># Finally, attributes are assigned to zones.  The first zone, &quot;default&quot;, con-</span>
<span class="c1"># tains most of the nodes; the zone &quot;bound&quot; contains only the inflow nodes.</span>
<span class="c1"># The initial water filling both zones is pure water, as signified by the</span>
<span class="c1"># asterisks in the water column, and the rock does not contain any species that</span>
<span class="c1"># participate in the reactions, as signified by the lack of a rock column.</span>
<span class="c1"># (The water column could also have been left out entirely, but is included</span>
<span class="c1"># here for clarity.)  No inflow occurs in the default zone, as shown by the</span>
<span class="c1"># asterisk in the boun column, but water of the  &quot;inflow&quot; type defined in the</span>
<span class="c1"># water block is flowing in through the nodes in the bound zone from time 1</span>
<span class="c1"># through time 4.167 of the simulation.  The dispersivity model &quot;dmod&quot; and the</span>
<span class="c1"># sorption model &quot;smod&quot; are applied to both zones.</span>
<span class="n">assign</span>       <span class="n">water</span>   <span class="n">boun</span>    <span class="n">time</span>    <span class="n">disp</span>    <span class="n">sorp</span>
<span class="n">default</span>      <span class="o">*</span>       <span class="o">*</span>       <span class="mi">0</span>       <span class="n">dmod</span>    <span class="n">smod</span>
<span class="n">bound</span>        <span class="o">*</span>       <span class="n">inflow</span>  <span class="mi">1</span><span class="o">&gt;</span><span class="mf">4.167</span> <span class="n">dmod</span>    <span class="n">smod</span>

<span class="n">end</span> <span class="n">trxn</span>
</pre></div>
</div>
<p>Below are the <code class="docutils literal notranslate"><span class="pre">trac</span></code> and <code class="docutils literal notranslate"><span class="pre">rxn</span></code> macros that were replaced by the above <code class="docutils literal notranslate"><span class="pre">trxn</span></code>, for reference and comparison.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span>trac
1.d-80 1.0 1.e-6 0.5
1. 2000 1.0 2000
5 5.0 1.e-6 2.8333e-3
6
1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 3.1623e-5 1.0 4.16667

1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 1.e-13 1.0 4.16667

1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 3.1623e-5 1.0 4.16667

0


0


0


end trac
rxn
** NCPLX, NUMRXN
2,4
** Coupling of the aqueous components (dRi/dUj)
2
1 0 1
0 1 0
** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
1    Cobalt[aq]     0     0     1.e-9
2      Iron[aq]     0     0     1.e-9
3      EDTA[aq]     0     0     1.e-9
** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
101   Co-EDTA[aq]   0
102   Fe-EDTA[aq]   0
** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immoblie spec, NIMM rows)
1   Co-EDTA[s]     0
2   Fe-EDTA[s]     0
3    Cobalt[s]     0
** IDVAP(IV), VAPNAM(IM), VAPPRT(IV) (ID # and name of vapor spec, NVAP rows)
** Skip nodes
0
** RSDMAX
1.0e-10
**** Chemical reaction information for equilibrium reactions ******
** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
         0
** CKEQ(IX) ,HEQ(IX) (Stability constants and Enthaplys, NCPLX rows)
1.0e+18    0
6.31e+27   0
** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
1.0       0.0       1.0
0.0       1.0       1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
      101      1
** Distribution coeffienct (kg water/ kg rock) **
    0.533
** Mass transfer coefficient (1/hr) **
      1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
        1      3
** Distribution coeffienct (kg rock/ kg water) **
     5.07
** Mass transfer coefficient (1/hr) **
      1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
      102      2
** Distribution coeffienct (kg rock/ kg water) **
    0.427
** Mass transfer coefficient (1/hr) **
      1.0
** GENERAL EXCHANGE REACTION (type 3) **
        3
** Where does the reaction take place? **
 1 0 0

** # of solid, liquid and vapor species **
        3   0   0
** forward and reverse rate constants (1/hr) **
  1.26e-2 0
** Solid Species in reaction **
  1      2      3
** Stoichiometry **
  1.0   -1.0  -1.0
end rxn
</pre></div>
</div>
</div>
<div class="section" id="additional-notes">
<h2>Additional Notes<a class="headerlink" href="#additional-notes" title="Permalink to this headline">¶</a></h2>
</div>
<div class="section" id="removed-features">
<h2>Removed Features<a class="headerlink" href="#removed-features" title="Permalink to this headline">¶</a></h2>
<p>Please note that the following features present in <cite>trac</cite> and <cite>rxn</cite> have been removed in <cite>trxn</cite>:</p>
<ul class="simple">
<li><code class="docutils literal notranslate"><span class="pre">trac</span></code></li>
</ul>
<blockquote>
<div><ul class="simple">
<li>Inflow concentrations for solids</li>
<li>Varying molecular diffusion coefficient by either species or location</li>
<li>Setting different dispersivities for different species</li>
</ul>
</div></blockquote>
<ul class="simple">
<li><code class="docutils literal notranslate"><span class="pre">rxn</span></code></li>
</ul>
<blockquote>
<div><ul class="simple">
<li>Using <cite>IFXCONC</cite> to specify that concentrations are free-ion only (all concentrations in <cite>water</cite>, <cite>rock</cite>, and <cite>gas</cite> must be total aqueous concentrations)</li>
<li>Enabling reactions at specific nodes only</li>
<li>Reaction type 6</li>
</ul>
</div></blockquote>
<ul class="simple">
<li>Both</li>
</ul>
<blockquote>
<div><ul class="simple">
<li>JA JB JC format for entering specific nodes for tracer injection, etc.  (These must be defined by zones.)</li>
</ul>
</div></blockquote>
<p>If these features are desired, old-style input using <cite>trac</cite> and <cite>rxn</cite> must be used instead.</p>
</div>
<div class="section" id="further-resources-and-verification">
<h2>Further resources and verification<a class="headerlink" href="#further-resources-and-verification" title="Permalink to this headline">¶</a></h2>
<p>Several test problems from the standard FEHM test suite that have <cite>trac</cite> and/or <cite>rxn</cite> macros have been converted to the <cite>trxn</cite> format.  These can be found in <cite>/scratch/nts/ms/trxn/test-problems</cite>.  The input files contain the original <cite>trac</cite> and/or <cite>rxn</cite> macros (turned off), along with an equivalent <cite>trxn</cite> macro.  The output for the <cite>trxn</cite> macro has been verified against the output for the original macros in each of these cases; the results of the comparison can be found in the <cite>plot/plot.png</cite> and <cite>plot/plot-orig.png</cite> files in the test directory.  The location of the input file in the test directory is included in parentheses after the test problem name.  The following test problems have been verified to work successfully with <cite>trxn</cite>:</p>
<ul class="simple">
<li><cite>baro_trans</cite> (<cite>input/baro_trans.in</cite>)</li>
<li><cite>cden_test (`input/static-multi1.dat</cite>)`</li>
<li><cite>dissolution</cite> (<cite>input/dissolution.in</cite>)</li>
<li><cite>fracture_transport</cite> (<cite>input/tangtestN.in</cite>)</li>
<li><cite>henrys_law</cite> (<cite>input/henrytest.in</cite>)</li>
<li><cite>multi_solute</cite> (<cite>input/multi_solute.in</cite>)</li>
<li><cite>sorption</cite> (<cite>input/sorption.in</cite>)</li>
<li><cite>transport3D</cite> (<cite>input/3d_trac.dat</cite>)</li>
</ul>
<p>In addition, several new test problems have been developed in order to better test the full range of <cite>trxn</cite>’s functionality.  These are also found in <cite>/scratch/nts/ms/trxn/test-problems</cite>, and documentation of the problem setup and expected results can be found in <cite>/scratch/nts/ms/trxn/test-problems/meta/doc</cite>.  The following test problems fall into this category:
* <cite>multirock</cite>
* <cite>inflow</cite>
* <cite>decay</cite></p>
</div>
<div class="section" id="bugs-error-handling-and-limitations">
<h2>Bugs, error handling, and limitations<a class="headerlink" href="#bugs-error-handling-and-limitations" title="Permalink to this headline">¶</a></h2>
<ul class="simple">
<li><cite>trxn</cite> will try to print out a nice error message if something goes wrong.  However, this is not guaranteed.</li>
<li>A <cite>zone</cite> macro with at least one valid zone is required before <cite>trxn</cite>.  This <cite>zone</cite> macro may contain no more than 100 zones.</li>
<li>If a block is specified multiple times, the values from the last block should be used; however, do not rely on this feature.</li>
<li>The maximum permitted length of input lines is 200 characters.</li>
<li>The maximum number of characters in a name (of a component, zone, model, etc.) is 40.</li>
<li>The maximum number of specifications (models, components, etc.) that any given block may contain is 100.</li>
<li>The zone macro preceding <cite>trxn</cite> may not contain more than 100 zones.</li>
<li>Model and species names should be strictly alphanumeric plus the five characters “(“, “)”, “+”, “_”, and “-“. Use of other characters may cause incorrect behavior.</li>
<li>If only the <cite>trac</cite> blocks of the macro are given, <cite>trxn</cite> should be compatible with an old <cite>rxn</cite> macro, as long as the <cite>trxn</cite> macro is read before the <cite>rxn</cite> macro.  However, this has not been tested and may not work reliably.  Furthermore, <cite>trxn</cite> is not compatible with the old <cite>trac</cite> macro.</li>
</ul>
</div>
<div class="section" id="debug-tools">
<h2>Debug tools<a class="headerlink" href="#debug-tools" title="Permalink to this headline">¶</a></h2>
<ul class="simple">
<li>The keyword <cite>debug</cite> in the <cite>ctrl</cite> block will enable some informational output that may be useful for debugging problems.  The <cite>stop</cite> keyword in the <cite>ctrl</cite> block will halt FEHM immediately after reading and processing <cite>trxn</cite>.</li>
<li><cite>null</cite> blocks are ignored by default, but are printed if debugging output is enabled.</li>
</ul>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ndar</span></code><a class="headerlink" href="#ndar" title="Permalink to this headline">¶</a></h1>

<p>The <span class="pre">ndar</span> macro enables non-Darcy flow.</p> 

<p>The use of Darcy’s law in subsurface porous flow models is nearly universal.  However, there are important applications where Darcy's law is not valid. 
  Flow in fractures and flow near pumping or extraction wells are applications where Non-Darcy flow models are often required. </p>

  

## Example macro ndar

The macro **ndar** is used to set non-darcy flow. The keyword *off* can be used to skip the macro and flow will be Darcy as usual. 

In this example, ```1 441 1``` is the node number start,stop,stride. ```1.0d9``` is the flow model which can be 1.d10, 1.d06, 1.0d9, and Beta default 1.d-15. Note for this test, 1.d06 is similar to the Darcy result, and 1.0d9 is slightly different (red lines in Figure 1).

<pre>
ndar
1 441 1 1.0d9

end ndar

or
ndar OFF (keyword **OFF** is used to skip a macro)
</pre>

## Example macro cont with vtk

The vtk option is available through the macro **cont**. In this example, vtk contour files are written at Time 0.0 and Time 50.0  *rootname*.00001_sca_node.vtk and *rootname*.00002_sca_node.vtk, each with contour data for Liquid Pressure (MPa) and Saturation. The pressure contour images were created using Time 0002 vtk files.

<pre>
cont
vtk     5000    1.00000e+19

</pre>


## Reference


Zeng, Z., Grigg, R. A Criterion for Non-Darcy Flow in Porous Media. Transp Porous Med 63, 57–69 (2006). [https://doi.org/10.1007/s11242-005-2720-3](https://doi.org/10.1007/s11242-005-2720-3)

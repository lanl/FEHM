try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

mat_prop = LegacyVTKReader( FileNames=[
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp_mat.vtk',
] )
RenameSource("model", mat_prop)
rv = GetRenderView()
dr = Show()
dr.ScalarOpacityUnitDistance = 1.7320508075688779
dr.EdgeColor = [0.0, 0.0, 0.5]

rv.CenterOfRotation = [   5.00000,    5.00000,    5.00000]

rv.CameraViewUp = [-0.4, -0.11, 0.92]
rv.CameraPosition = [  30.00000,   20.00000,   20.00000]
rv.CameraFocalPoint = [   5.00000,    5.00000,    5.00000]

mr = GetDisplayProperties(mat_prop)
mr.Representation = 'Surface With Edges'

lt = GetLookupTableForArray( 'perm_x', 1, RGBPoints=[-20.00, 0.23, 0.299, 0.754, -14.00, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

pf = CreatePiecewiseFunction( Points=[-20.00, 0.0, 0.5, 0.0, -14.00, 1.0, 0.5, 0.0] )

mr.ScalarOpacityFunction = pf
mr.ColorArrayName = ('POINT_DATA', 'perm_x')
mr.LookupTable = lt

lt.ScalarOpacityFunction = pf

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='perm_x', LabelFontSize=12, Enabled=1, TitleFontSize=12 )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

lt = GetLookupTableForArray('perm_x', 1 )

ScalarBarWidgetRepresentation1.LookupTable = lt

AnimationScene1 = GetAnimationScene()
AnimationScene1.AnimationTime = 0.0
rv.ViewTime = 0.0
source = FindSource("model")
SetActiveSource(source)

G = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )
G.GlyphTransform = "Transform2"
G.GlyphType = "Sphere"
G.RandomMode = 0
G.ScaleMode = 'off'
G.MaskPoints = 0
G.GlyphType.Radius =    0.10000

RenameSource("nodes", G)

rv = GetRenderView()
mr = GetDisplayProperties(source)
dr = Show()
dr.ColorArrayName = ('POINT_DATA', 'n')
dr.ScaleFactor = 1.1
dr.SelectionPointFieldDataArrayName = "nodes"
dr.EdgeColor = [0.0, 0.0, 0.5000076295109483]
dr.ColorArrayName = ('POINT_DATA', '')
dr.DiffuseColor = [0.,0.,0.]
dr.Visibility = 0
cols = [
[1.00,1.00,0.00],
[1.00,0.00,1.00],
[0.00,1.00,1.00],
]
zones = [
'zone0001_lower',
'zone0002_middle',
'zone0003_upper',
]
for zone,col in zip(zones,cols):
	AnimationScene1 = GetAnimationScene()
	AnimationScene1.AnimationTime = 0.0
	rv.ViewTime = 0.0
	source = FindSource("model")
	SetActiveSource(source)
	
	G = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )
	G.GlyphTransform = "Transform2"
	G.Scalars = ['POINTS', zone]
	G.ScaleMode = 'scalar'
	G.GlyphType = "Sphere"
	G.RandomMode = 0
	G.MaskPoints = 0
	
	G.GlyphType.Radius =    0.20000
	
	RenameSource(zone, G)
	
	rv = GetRenderView()
	mr = GetDisplayProperties(source)
	dr = Show()
	dr.ColorArrayName = ('POINT_DATA', 'n')
	dr.ScaleFactor = 1.1
	dr.SelectionPointFieldDataArrayName = zone
	dr.EdgeColor = [0.0, 0.0, 0.5000076295109483]
	dr.Opacity = 0.5
	
	lt = GetLookupTableForArray(zone, 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.5, 0.865, 0.865, 0.865, 1.0]+col, VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )
	
	pf = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
	
	dr.ColorArrayName = ('POINT_DATA', zone)
	dr.LookupTable = lt
	dr.Visibility = 0
	
	lt.ScalarOpacityFunction = pf
contour_output = LegacyVTKReader( FileNames=[
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0000.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0001.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0002.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0003.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0004.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0005.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0006.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0007.vtk',
'/scratch/er/dharp/source/pyfehm/tutorials/tut1/temp.0008.vtk',
] )
RenameSource("contour_output", contour_output)
lt = GetLookupTableForArray('P', 1, RGBPoints=[   4.00000, 0.23, 0.299, 0.754,    6.11555, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

pf = CreatePiecewiseFunction( Points=[   4.00000, 0.0, 0.5, 0.0,    6.11555, 1.0, 0.5, 0.0] )

dr = Show() #dr = DataRepresentation1
dr.Representation = 'Surface With Edges'
dr.EdgeColor = [0.15, 0.15, 0.15]
dr.ScalarOpacityFunction = pf
dr.ColorArrayName = ('POINT_DATA', 'P')
dr.ScalarOpacityUnitDistance = 1.7320508075688779
dr.LookupTable = lt

rv.ViewTime =    9

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='P', LabelFontSize=12, Enabled=1, LookupTable=lt, TitleFontSize=12 )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

model = FindSource("model")
model_rep = GetDisplayProperties(model)
contour_output = FindSource("contour_output")
cont_rep = GetDisplayProperties(contour_output)
model_rep.Visibility = 1
cont_rep.Visibility = 0	
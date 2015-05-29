import vtk
import numpy as num 
from pyrocko import model
from pyrocko import orthodrome as ortho
import os

def make_polydata_actor(centers, normals, colors, return_pdm=False):
    """ Create the actor and set colors

    :param return_pdm: if True give back the polydatamapper
    :param centers: list of centers as tuples 
    :param normals: list of normals as tuples
    :param colors: list of rgb tuples
    """
     
    # In order to build the bounding box it is convenient to combine
    # All sources using a vtkAppenPolyData
    apd = vtk.vtkAppendPolyData()

    mappers = []
    # create source
    for i in range(len(centers)):

        normal = normals[i]
        if normal is None:
            continue
        source = vtk.vtkRegularPolygonSource()
        source.SetNormal(normal)
        source.SetCenter(*centers[i])
        source.SetNumberOfSides(16)
        source.SetRadius(100)
        apd.AddInput(source.GetOutput())
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(apd.GetOutputPort())
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    # set colors:
    actor.GetProperty().SetColor(*colors[i])
    if return_pdm:
        return actor, apd
    else:
        return actor

    
def setup_renderer(renderer, actors, bboxpolydata=None):
    ren = renderer
    for actor in actors:
        # assign actor to the renderer
        ren.AddActor(actor)
    ############################################
    # Create a vtkOutlineFilter to draw the bounding box of the data set.
    # Also create the associated mapper and actor.
    if bboxpolydata is not None:
        outline = vtk.vtkOutlineFilter()
        outline.SetInput(bboxpolydata.GetOutput())
        mapOutline = vtk.vtkPolyDataMapper()
        mapOutline.SetInputConnection(outline.GetOutputPort())
        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(mapOutline)
        outlineActor.GetProperty().SetColor(100, 100, 100)
        ren.AddViewProp(outlineActor)
        ren.AddActor(outlineActor)
    

    # Create a vtkLight, and set the light parameters.
    #light = vtk.vtkLight()
    #light.SetFocalPoint(0.21406, 1.5, 0)
    #light.SetPosition(8.3761, 4.94858, 4.12505)
    #ren.AddLight(light)
    return ren

def moment_tensors2normals(tensors):
    """ NEEDS REVISION and MODIFICATION"""
    normals = []
    for mt in tensors:
        if mt:
            normals.append(mt.p_axis().tolist()[0])
        else:
            normals.append(None)

    return normals

def to_cartesian(items):
    res = []
    latlon00 = ortho.Loc(0.,0.)
    for i, item in enumerate(items):

        y, x = ortho.latlon_to_ne(latlon00, item)
        depth = item.depth
        lat = item.lat/180.*num.pi
        res.append((x, y, -depth))
    return res

def to_colors(items):
    """2 implement"""
    return [(1,0,0)]*len(items)

def read_data(event_fn=None, events=None):
    if event_fn is not None and events is None:
        events = model.load_events(event_fn)
    else:
        events = events
    moment_tensors = []
    for e in events:
        if e.moment_tensor is not None:
            moment_tensors.append(e.moment_tensor)
        else:
            moment_tensors.append(None)

    normals = moment_tensors2normals(moment_tensors) 
    centers = to_cartesian(events)
    colors = to_colors(moment_tensors)
    
    return normals, centers, colors

if __name__=="__main__":
    webnet = os.environ['WEBNET']
    input_fn = webnet+"/meta/events2008_mt.pf"
    compare_fn = webnet+"/rapid_compile/rapidinv_events.pf"
    ren = vtk.vtkRenderer()
    actors = [] 
    normals, centers, colors = read_data(compare_fn)
    colors = [(0,1,0)]*len(colors)
    kwargs = {"centers": centers, 'normals':normals, "colors":colors, "return_pdm":True}
    actor1, apd = make_polydata_actor(**kwargs)
    actors.append(actor1)

    #normals, centers, colors = read_data(input_fn)
    #kwargs = {"centers": centers, 'normals':normals, "colors":colors, "return_pdm":True}
    #actor2, apd = make_polydata_actor(**kwargs)
    #actors.append(actor2)
    #
    
    #setup_renderer(ren,[actor1, actor2], bboxpolydata=apd)
    setup_renderer(ren, actors, bboxpolydata=apd)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
     
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
         
    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

    

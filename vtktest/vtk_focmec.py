import vtk
import numpy as num 
from pyrocko import model, moment_tensor
from pyrocko import orthodrome as ortho
import os


def get_fault_planes(p_axis, t_axis, null_axis):
    pplanes = []
    tplanes = []
    for i in range(len(null_axis)):
        rotmat = moment_tensor.rotation_from_axis_and_angle(angle=45, axis=null_axis[i])
        tplanes.append(num.array(t_axis[i]*rotmat)[0])
        pplanes.append(num.array(p_axis[i]*rotmat)[0])
    return [pplanes, tplanes]

def make_polydata_actor(centers, normals, return_pdm=False, type='circle'):
    """ Create the actor and set colors

    :param return_pdm: if True give back the polydatamapper
    :param centers: list of centers as tuples 
    :param normals: list of list of normals as tuples
    """
     
    # In order to build the bounding box it is convenient to combine
    # All sources using a vtkAppenPolyData
    apd = vtk.vtkAppendPolyData()
    #save_coordinates(centers)
    mean_center, plane_normal = plane_fit(centers)
    source = vtk.vtkRegularPolygonSource()
    source.SetNumberOfSides(16)
    source.SetRadius(800)
    source.SetNormal(*plane_normal)
    apd.AddInput(source.GetOutput())

    mappers = []
    # create source
    for i in range(len(centers.T)):
        normal = normals[i]
        if normal is None:
            continue
        if type=='torus':
            source = vtk.vtkSuperquadricSource();
            source.SetScale(1.0, 1.0, 1.0)
            source.SetPhiResolution (16)
            source.SetThetaResolution(16)
            source.SetThetaRoundness (1)
            source.SetThickness (0.1)
            source.SetSize(100.5)
            source.SetToroidal(1) 
            source.SetNormal(normal)
        elif type=='circle':
            source = vtk.vtkRegularPolygonSource()
            source.SetNumberOfSides(16)
            #source.SetRadius(100)
            source.SetRadius(10)
            source.SetNormal(normal)
        source.SetCenter(*centers.T[i])
        apd.AddInput(source.GetOutput())
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(apd.GetOutputPort())

    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
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

def moment_tensors2normals(tensors, get):
    """ NEEDS REVISION and MODIFICATION"""
    normals = []
    for mt in tensors:
        if mt:
            n = getattr(mt, get)()
            normals.append(n.tolist()[0])
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
    res = num.array(res) 
    res = res.T
    return res

def to_colors(items):
    """2 implement"""
    return [(1,0,0)]*len(items)

def read_data(event_fn=None, events=None, get=None):
    
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
    
    normals = []
    for g in get:
        normals.append(moment_tensors2normals(moment_tensors, g))
    centers = to_cartesian(events)
    colors = to_colors(moment_tensors)
    
    return normals, centers, colors

def save_coordinates(coordinates):
    savestr = ''
    for x,y,z in coordinates.T:
        savestr += '%s  %s  %s\n' % (x,y,z)

    with open('coordinates.dat', 'w') as f:
        f.write(savestr)

def plane_fit(coordinates):
    '''Calculate a best fitting plane through the point cloud. 
    Returns the normal and the mean center.'''
    data = coordinates
    cm = num.mean(data, 1)

    # de-mean data
    data[0] -= cm[0]
    data[1] -= cm[1]
    data[2] -= cm[2]

    u, s, vh = num.linalg.svd(data.T)
    vh = vh.conj().transpose()
    return cm, vh[:, -1]

def make_plane(normal):
    plane = vtk.vtkPlaneSource()
    plane.SetNormal(normal)
    return plane

if __name__=="__main__":
    webnet = os.environ['WEBNET']
    input_fn = webnet+"/meta/events2008_mt.pf"
    compare_fn = webnet+"/rapid_compile/rapidinv_events.pf"
    ren = vtk.vtkRenderer()
    actors = [] 
    normals_list, centers, colors = read_data(compare_fn, get=['p_axis', 't_axis', 'null_axis'])
    plot_fault_planes = True

    if plot_fault_planes:
        normals = get_fault_planes(*normals_list)
    else:
        normals = normals_list[:2]
    
    for i in range(2):
        kwargs = {"centers": centers, 'normals': normals[i], "return_pdm":True}
        if i==0:
            color = (0,1,0)
            opacity = (1.)
        else:
            color = (1,1,1)
            opacity = (0.1)
    
        actor1, apd = make_polydata_actor(**kwargs)
        actor1.GetProperty().SetColor(color)
        actor1.GetProperty().SetOpacity(opacity)
        actors.append(actor1)

    #normals_list, centers, colors = read_data(input_fn, get=['p_axis', 't_axis', 'null_axis'])
    #if plot_fault_planes:
    #    normals = get_fault_planes(*normals_list)
    #else:
    #    normals = normals_list[:2]
     
    #for i in range(2):
    #    kwargs = {"centers": centers, 'normals':normals[i], "return_pdm":True}
    #    if i==0:
    #        color = (1,0,0)
    #        opacity = (1.)
    #    else:
    #        color = (1,1,1)
    #        opacity = (0.1)
    #    actor2, apd = make_polydata_actor(**kwargs)
    #    actor2.GetProperty().SetColor(color)
    #    actor2.GetProperty().SetOpacity(opacity)
    #    actors.append(actor2)
    #
    #
    setup_renderer(ren, actors, bboxpolydata=apd)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
     
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
         
    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

    

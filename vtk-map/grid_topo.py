#!/usr/bin/env python
import time

import vtk
import numpy as num
from pyrocko import gmtpy
from pyrocko.automap import Map
from pyrocko import orthodrome as od


def get_topography_actor(x, y, topography, super_elevation, decimation,
                         smoothing=False):
    topography = topography.T
    topography *= super_elevation

    points, triangles, colors = setup_topography(x, y, topography,
                                                 decimation=decimation)
    # Create a polydata object
    trianglePolyData = vtk.vtkPolyData()

    # Add the geometry and topology to the polydata
    trianglePolyData.SetPoints(points)
    trianglePolyData.GetPointData().SetScalars(colors)
    trianglePolyData.SetPolys(triangles)

    # Clean the polydata so that the edges are shared !
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInput(trianglePolyData)

    # Use a filter to smooth the data (will add triangles and smooth)
    # Use two different filters to show the difference
    mapper = vtk.vtkPolyDataMapper()

    if smoothing:
        smooth_loop = vtk.vtkLoopSubdivisionFilter()
        smooth_loop.SetNumberOfSubdivisions(2)
        smooth_loop.SetInputConnection(cleanPolyData.GetOutputPort())
        mapper.SetInputConnection(smooth_loop.GetOutputPort())

    else:
        mapper.SetInput(trianglePolyData)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


def scaled_mesh(xin, yin):
    x, y = num.meshgrid(xin, yin)
    xs = num.zeros(x.shape)
    ys = num.zeros(y.shape)
    olats = num.zeros(y.shape[1])
    olons = num.zeros(x.shape[1])

    for i in xrange(len(x)):
        xs[i], ys[i] = od.latlon_to_ne_numpy(olats, olons, y[i], x[i])

    return xs, ys


def setup_topography(x, y, topography, xmax=None, ymax=None, decimation=1):
    # Define points, triangles and colors
    x = x[::decimation]
    y = y[::decimation]
    lonsize = len(x)-1 if not xmax else xmax
    latsize = len(y)-1 if not ymax else ymax
    colors = vtk.vtkUnsignedCharArray()
    # colors.SetNumberOfComponents(3)
    colors.SetNumberOfComponents(1)
    points = vtk.vtkPoints()
    triangles = vtk.vtkCellArray()
    zmax = topography.max()
    zmin = topography.min()
    zrange = zmax - zmin
    xmesh, ymesh = scaled_mesh(x, y)
    count = 0
    t1 = time.time()
    topography = topography.T
    topo_new = num.zeros((len(y), len(x)))
    for iy in xrange(len(y)):
        topo_new[iy, :] = topography[iy*decimation, ::decimation]
    topography = topo_new
    for i in xrange(latsize):
        print '%i / %i' % (i+1, latsize)
        for j in xrange(lonsize-3):

            d = (ymesh[i][j], xmesh[i][j], topography[i][j])
            c = (ymesh[i][j+1], xmesh[i][j+1], topography[i][j+1])
            b = (ymesh[i+1][j+1], xmesh[i+1][j+1], topography[i+1][j+1])
            a = (ymesh[i+1][j], xmesh[i+1][j], topography[i+1][j])
            points.InsertNextPoint(*a)
            points.InsertNextPoint(*b)
            points.InsertNextPoint(*c)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count)
            triangle.GetPointIds().SetId(1, count + 1)
            triangle.GetPointIds().SetId(2, count + 2)

            triangles.InsertNextCell(triangle)

            points.InsertNextPoint(*a)
            points.InsertNextPoint(*d)
            points.InsertNextPoint(*c)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, count + 3)
            triangle.GetPointIds().SetId(1, count + 4)
            triangle.GetPointIds().SetId(2, count + 5)

            count += 6

            triangles.InsertNextCell(triangle)

            # rs = [[int((zmax-topography[j][i])/zrange*255)]]*6
            rs = [[int((zmax-topography[i][j])/zrange*255)]]*6
            map(colors.InsertNextTupleValue, rs)
    print 'total time needed ', time.time()-t1
    return points, triangles, colors


def setup_vtk_map_actor(lat, lon, radius, super_elevation=1, decimation=1,
                        smoothing=False):
    m = Map(lat=lat, lon=lon, radius=radius)
    m._setup()
    fn, ilum = m._prep_topo('land')
    topodata = gmtpy.loadgrd(fn)
    return get_topography_actor(*topodata, super_elevation=super_elevation,
                                decimation=decimation, smoothing=smoothing)


if __name__ == '__main__':

    lat = 46.
    lon = 10.
    radius = 25000
    actor = setup_vtk_map_actor(lat, lon, radius, decimation=5, smoothing=True)
    #
    # mapper = vtk.vtkPolyDataMapper()
    # mapper.SetInputConnection(smooth_butterfly.GetOutputPort())
    # actor_butterfly = vtk.vtkActor()
    # actor_butterfly.SetMapper(mapper)
    # actor_butterfly.SetPosition(64, 0, 0)

    # Visualise
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Add actors and render
    renderer.AddActor(actor)

    renderer.SetBackground(1, 1, 1)  # Background color white
    renderWindow.SetSize(800, 800)
    renderWindow.Render()
    renderWindowInteractor.Start()

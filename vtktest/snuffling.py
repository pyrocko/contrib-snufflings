import numpy as num
import vtk
from vtk.util import numpy_support

from pyrocko.snuffling import Snuffling, Param, Switch
from pyrocko import orthodrome as ortho

class VtkTest(Snuffling):

    ''' Create 3D interactive plots of events and stations using vtk.

    This snuffling requires VTK to be installed including Python wrappers.
    On debian using apt-get:
    apt-get install python-vtk
    '''

    def setup(self):
        self.set_name('VTK')
        self.add_parameter(Param('Superelevation', 'z_scale', 1., 1., 1000.))
        self.add_parameter(Switch('Stations', 'want_stations', True))
        self.add_parameter(Switch('Events', 'want_events', True))
        self.actors = []
        self.frame = None
        self.scale_sources = 1.

    def locations_to_ned(self, locations):
        npoints = len(locations)
        lats = num.zeros(npoints)
        lons = num.zeros(npoints)
        depths = num.zeros(npoints)

        for i_e, e in enumerate(locations):
            lats[i_e] = float(e.lat)
            lons[i_e] = float(e.lon)
            depths[i_e] = float(e.depth)

        depths *= self.z_scale
        nz = num.zeros(npoints)
        ns, es = ortho.latlon_to_ne_numpy(nz, nz, lats, lons)

        return ns, es, depths

    def stations_to_vtkcone_actors(self, data):
        actors = []
        for i in xrange(data.GetNumberOfTuples()):
            source = vtk.vtkConeSource()
            source.SetCenter(*data.GetTuple3(i))
            source.SetRadius(200)
            source.SetHeight(200)
            source.Update()
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actors.append(actor)

        return actors

    def events_to_vtksphere_actors(self, data):
        points = vtk.vtkPoints()
        actors = []
        npoints =None
        for i in xrange(data.GetNumberOfTuples()):
            source = vtk.vtkSphereSource()
            source.SetCenter(*data.GetTuple3(i))
            source.SetRadius(50)
            source.Update()
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actors.append(actor)

        return actors

    def call(self):
        self.cleanup()
        viewer = self.get_viewer()
        if self.want_events:
            markers = self.get_selected_event_markers()
            events = [m.get_event() for m in markers]

            ns, es, depths = self.locations_to_ned(events)
            adata = num.array((depths, ns, es))
            valrange_n = (num.min(adata[1, :]), num.max(adata[1, :]))
            valrange_e = (num.min(adata[2, :]), num.max(adata[2, :]))

            adata = adata.flatten(order='F')
            data = numpy_support.numpy_to_vtk(
                adata, deep=True, array_type=vtk.VTK_FLOAT )
            data.SetNumberOfComponents(3)

            sphere_actors = self.events_to_vtksphere_actors(data)
        else:
            sphere_actors = []

        if self.want_stations:
            stations = self.get_stations()

            ns, es, depths = self.locations_to_ned(stations)
            adata = num.array((depths, ns, es)).flatten(order='F')
            data = numpy_support.numpy_to_vtk(
                adata, deep=True, array_type=vtk.VTK_FLOAT)
            data.SetNumberOfComponents(3)

            cone_actors = self.stations_to_vtkcone_actors(data)
        else:
            cone_actors = []

        frame = self.vtk_frame()

        for actor in cone_actors:
            frame.add_actor(actor)

        for actor in sphere_actors:
            frame.add_actor(actor)
        frame.renderer.SetBackground(0.1, 0.2, 0.4)
        frame.init()


def __snufflings__():
    return [VtkTest()]


import numpy as num

from pyrocko.snuffling import Snuffling, Param, Switch
from pyrocko import orthodrome as ortho


class VtkTest(Snuffling):

    ''' Create 3D interactive plots of events and stations

    This snuffling requires VTK to be installed including Python wrappers.
    On debian using apt-get:
    apt-get install python-vtk

    You need access to the SRTMGL3 database from

    https://earthdata.nasa.gov

    which requires username and password. If you do not have an account, yet
    follow the given link and create one. Username and password have
    to be added to be your ~/.pyrocko/config.pf as follows:

    earthdata_credentials: [username, password]
    '''

    def setup(self):
        self.set_name('VTK 3D-Map')
        self.add_parameter(Param(
            'Vertical exaggeration', 'z_scale', 1., 1., 1000.))
        self.add_parameter(Param(
            'Topographic decimation', 'z_decimation', 1, 1, 12,
            low_is_none=True))
        self.add_parameter(Param(
            'Margin Radius [km]', 'margin_radius', 1, 1, 100))
        self.add_parameter(Switch('Stations', 'want_stations', True))
        self.add_parameter(Switch('Events', 'want_events', True))
        self.add_parameter(Switch('Topography', 'want_topo', True))
        self.add_parameter(Switch('Topography smoothing', 'smoothing', False))
        self.actors = []
        self.topo_actor = None
        self.frame = None
        self.set_live_update(False)

    def locations_to_ned(self, locations, has_elevation=False):
        npoints = len(locations)
        lats = num.zeros(npoints)
        lons = num.zeros(npoints)
        depths = num.zeros(npoints)

        for i_e, e in enumerate(locations):
            lats[i_e] = float(e.lat)
            lons[i_e] = float(e.lon)
            if has_elevation:
                depths[i_e] = float(e.depth) - float(e.elevation)
            else:
                depths[i_e] = float(e.depth)

        depths *= self.z_scale
        nz = num.zeros(npoints)
        ns, es = ortho.latlon_to_ne_numpy(nz, nz, lats, lons)

        return ns, es, depths

    def stations_to_vtkcone_actors(self, data, size=500.):
        actors = []
        for i in xrange(data.GetNumberOfTuples()):
            source = vtk.vtkConeSource()
            source.SetCenter(*data.GetTuple3(i))
            source.SetRadius(size)
            source.SetHeight(size)
            source.SetDirection(0., 0., -1.)
            source.Update()
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actors.append(actor)

        return actors

    def events_to_vtksphere_actors(self, data, size=500.):
        actors = []
        for i in xrange(data.GetNumberOfTuples()):
            source = vtk.vtkSphereSource()
            source.SetCenter(*data.GetTuple3(i))
            source.SetRadius(size)
            source.Update()
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actors.append(actor)

        return actors

    def call(self):
        try:
            global vtk
            import vtk
            from vtk.util import numpy_support
            import sys
            sys.path[0:0] = [self.module_dir()]
            from grid_topo import setup_vtk_map_actor
            sys.path[0:1] = []
        except ImportError as _import_error:
            self.fail('\nImportError:\n%s' % _import_error)
            vtk = None
        self.cleanup()
        stations = []
        events = []
        cone_actors = []
        sphere_actors = []

        if self.want_stations:
            stations = self.get_stations()

        if self.want_events:
            markers = self.get_selected_event_markers()
            events = [m.get_event() for m in markers]

        all_lats = []
        all_lons = []

        for s in stations:
            all_lats.append(s.lat)
            all_lons.append(s.lon)

        for e in events:
            all_lats.append(e.lat)
            all_lons.append(e.lon)

        center_lat, center_lon = ortho.geographic_midpoint(
            num.array(all_lats), num.array(all_lons))

        center_lats = num.array([center_lat]*len(all_lats))
        center_lons = num.array([center_lon]*len(all_lons))
        distances = ortho.distance_accurate50m_numpy(
            num.array(all_lats), num.array(all_lons),
            center_lats, center_lons)
        distance_max = num.max(distances)
        size = distance_max / 50.

        if len(events) != 0:
            ns, es, depths = self.locations_to_ned(events)
            adata = num.array((es, ns, -depths))
            adata = adata.flatten(order='F')
            data = numpy_support.numpy_to_vtk(
                adata, deep=True, array_type=vtk.VTK_FLOAT)
            data.SetNumberOfComponents(3)

            sphere_actors = self.events_to_vtksphere_actors(data, size=size/2.)

        if len(stations) != 0:
            ns, es, depths = self.locations_to_ned(stations,
                                                   has_elevation=True)
            adata = num.array((es, ns, -depths)).flatten(order='F')
            data = numpy_support.numpy_to_vtk(
                adata, deep=True, array_type=vtk.VTK_FLOAT)
            data.SetNumberOfComponents(3)

            cone_actors = self.stations_to_vtkcone_actors(data, size=size)

        if self.want_topo:
            distance_max += self.margin_radius * 1000
            self.topo_actor = setup_vtk_map_actor(
                center_lat, center_lon, distance_max,
                super_elevation=self.z_scale,
                decimation=int(self.z_decimation or 1),
                smoothing=self.smoothing)

        frame = self.vtk_frame()

        for actor in cone_actors:
            actor.GetProperty().SetColor(0., 0., 1.)
            frame.add_actor(actor)

        for actor in sphere_actors:
            actor.GetProperty().SetColor(1., 0., 0.)
            frame.add_actor(actor)

        if self.topo_actor:
            frame.add_actor(self.topo_actor)
        frame.renderer.SetBackground(0.1, 0.2, 0.4)
        frame.init()


def __snufflings__():
    return [VtkTest()]

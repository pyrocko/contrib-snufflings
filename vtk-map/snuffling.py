import numpy as num
import matplotlib
import matplotlib.pyplot as plt

from pyrocko.gui.snuffling import Snuffling, Param, Switch
from pyrocko import orthodrome as ortho

try:
    from vtk.util import numpy_support
except ImportError as e:
    vtk = None


class ColorMapper():
    def __init__(self, cmap):
        self.cmap = plt.get_cmap(cmap)
        self.set_range(0., 1.)

    def set_range(self, min, max):
        self.mapping = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(min, max),
            cmap=self.cmap)

    def __call__(self, v):
        return self.mapping.to_rgba(v)


def locations_to_ned(locations, has_elevation=False, z_scale=1.):
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

    depths *= z_scale
    nz = num.zeros(npoints)
    ns, es = ortho.latlon_to_ne_numpy(nz, nz, lats, lons)

    return ns, es, depths


def events_to_vtksphere_actors(events, z_scale, size=500.):

    ns, es, depths = locations_to_ned(events)
    times = [e.time for e in events]

    if not times:
        times = num.ones(ntuples)

    adata = num.array((es, ns, -depths))
    adata = adata.flatten(order='F')
    data = numpy_support.numpy_to_vtk(
        adata, deep=True, array_type=vtk.VTK_FLOAT)
    data.SetNumberOfComponents(3)
    actors = []
    ntuples = data.GetNumberOfTuples()

    cmap = plt.get_cmap

    for i in range(ntuples):
        source = vtk.vtkSphereSource()
        source.SetCenter(*data.GetTuple3(i))
        source.SetRadius(size)
        source.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        r, g, b, a = events[i].get_vtk_color(events[i].time)
        actor.GetProperty().SetColor(r, g, b)
        actors.append(actor)

    return actors


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
        self.add_trigger('Make screenshot', self.save_image)
        self.actors = []
        self.frame = None
        self.set_live_update(False)


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

    def call(self):
        try:
            global vtk
            import vtk
            import sys
            sys.path[0:0] = [self.module_dir()]
            from grid_topo import setup_vtk_map_actor
            sys.path[0:1] = []
        except ImportError as _import_error:
            self.fail('\nImportError:\n%s' % _import_error)
            vtk = None

        self.cleanup()
        viewer = self.get_viewer()

        stations = []
        events = []
        cone_actors = []
        sphere_actors = []

        if self.want_stations:
            stations = self.get_stations()

        if self.want_events:
            markers = self.get_selected_event_markers()
            if len(markers) == 0:
                tmin, tmax = self.get_selected_time_range(fallback=True)
                markers = filter(lambda x: tmin < x.tmin < tmax,
                                 self.get_event_markers())
            events = [m.get_event() for m in markers]

            active_event = viewer.get_active_event()
            to_rgba = ColorMapper('summer')
            times = [e.time for e in events]
            to_rgba.set_range(min(times), max(times))
            if active_event:


                to_rgba = ColorMapper('gray')
                to_rgba.set_range(min(times), max(times))

            for e in events:
                e.get_vtk_color = to_rgba

            if active_event:
                def return_red(e):
                    return (1., 0., 0., 1.)

                active_event.get_vtk_color = return_red
                events.append(active_event)

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
            sphere_actors = events_to_vtksphere_actors(
                events, z_scale=self.z_scale, size=size/2.)

        if len(stations) != 0:
            ns, es, depths = locations_to_ned(
                stations, z_scale=self.z_scale, has_elevation=True)
            adata = num.array((es, ns, -depths)).flatten(order='F')
            data = numpy_support.numpy_to_vtk(
                adata, deep=True, array_type=vtk.VTK_FLOAT)
            data.SetNumberOfComponents(3)

            cone_actors = self.stations_to_vtkcone_actors(data, size=size)

        if self.want_topo:
            distance_max += self.margin_radius * 1000.
            topo_actor = setup_vtk_map_actor(
                center_lat, center_lon, distance_max,
                super_elevation=self.z_scale,
                decimation=int(self.z_decimation or 1),
                smoothing=self.smoothing)

        self.frame = self.vtk_frame()

        for actor in cone_actors:
            actor.GetProperty().SetColor(0., 0., 1.)
            self.frame.add_actor(actor)

        for actor in sphere_actors:
            self.frame.add_actor(actor)

        if self.want_topo:
            self.frame.add_actor(topo_actor)
        self.frame.renderer.SetBackground(0.01, 0.05, 0.1)

        self.frame.init()

    def save_image(self):
        fn = self.output_filename('save PNG')
        renWin = self.frame.vtk_widget.GetRenderWindow()

        imageFilter = vtk.vtkWindowToImageFilter()
        imageFilter.SetInput(renWin)
        imageFilter.SetInputBufferTypeToRGB()
        imageFilter.ReadFrontBufferOff()
        imageFilter.Update()

        pngwriter = vtk.vtkPNGWriter()
        pngwriter.SetInputConnection(imageFilter.GetOutputPort())
        pngwriter.SetFileName(fn)
        pngwriter.Write()


def __snufflings__():
    return [VtkTest()]

import numpy as num

from PyQt4.QtCore import *
from PyQt4.QtGui import QFileDialog
from pyrocko import moment_tensor, model
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
from pyrocko import gf

km = 1000.

class Seismosizer(Snuffling):

    def __init__(self):
        Snuffling.__init__(self)

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Seismosizer')
        self.add_parameter(Param('Time', 'time', 0.0, -50., 50.))
        self.add_parameter(Param('Latitude', 'lat', 0.0, -90., 90.))
        self.add_parameter(Param('Longitude', 'lon', 0.0, -180., 180.))
        self.add_parameter(Param('North shift', 'north_km', 0.0, -50., 50.))
        self.add_parameter(Param('East shift', 'east_km', 0.0, -50., 50.))
        self.add_parameter(Param('Depth', 'depth_km', 10.0, 0.0, 600.0))
        self.add_parameter(Param('Magnitude', 'magnitude', 6.0, 0.0, 10.0))
        self.add_parameter(Param('Strike', 'strike', 0., -180., 180.))
        self.add_parameter(Param('Dip', 'dip', 90., 0., 90.))
        self.add_parameter(Param('Rake', 'rake', 0., -180., 180.))
        self.add_parameter(Param('Length', 'length', 0., 0., 1000*km))
        self.add_parameter(Param('Width', 'width', 0., 0., 500*km))
        self.add_parameter(Param('Nucleation X', 'nucleation_x', -1., -1., 1.))
        self.add_parameter(Param('Rise-time', 'risetime', 0.0, 0.0, 20.0))
        self.add_parameter(Choice('GF Store', 'store_id', '<not loaded yet>', ['<not loaded yet>']))
        
        self.add_trigger('Set Engine', self.set_engine)
        self.add_trigger('Set Params from Event', self.mechanism_from_event)
        self.add_trigger('Add Stores', self.add_store)

        self.store_ids = None
        self.offline_config = None
        self._engine = None

    def setup_gui(self, *args, **kwargs):
        def visibility(visible):
            if visible:
                if self._engine is None:
                    self.set_engine()

        retval = Snuffling.setup_gui(self, *args, **kwargs)
        panel = self._panel.parent()
        panel.connect( panel, SIGNAL('visibilityChanged(bool)'), visibility)

        return retval

    def set_engine(self):
        self._engine = None
        self.store_ids = self.get_store_ids()
        self.set_parameter_choices('store_id', self.store_ids)
        if 'global_2s' in self.store_ids:
            self.store_id = 'global_2s'
        else:
            self.store_id = self.store_ids[0]


    def get_engine(self):
        if not self._engine:
            self._engine = gf.LocalEngine(use_config=True)

        return self._engine

    def get_store_ids(self):
        return self.get_engine().get_store_ids()
        
    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()

        # get time range visible in viewer
        viewer = self.get_viewer()

        event = viewer.get_active_event()
        if event:
            event, stations = self.get_active_event_and_stations(missing='raise')
        else:
            event = model.Event(lat=self.lat, lon=self.lon)
            stations = []
    
        if not stations:
            stations = []
            for (lat, lon) in [(5.,0.), (-5.,0.)]:
                stations.append(
                    model.Station(station='(%g, %g)' % (lat, lon), lat=lat, lon=lon))
        
        source = gf.RectangularSource(
            time=event.time+self.time,
            lat=event.lat,
            lon=event.lon,
            north_shift=self.north_km*km,
            east_shift=self.east_km*km,
            depth=self.depth_km*km,
            magnitude=self.magnitude,
            strike=self.strike,
            dip=self.dip,
            rake=self.rake,
            length=self.length,
            width=self.width,
            nucleation_x=self.nucleation_x,
            risetime=self.risetime)

        source.regularize()

        targets = []
        for station in stations:
            for component in 'ZNE':

                target = gf.Target(
                    codes=(
                        station.network,
                        station.station,
                        '',
                        component),
                    lat=station.lat,
                    lon=station.lon,
                    store_id=self.store_id,
                    optimization='enable',
                    interpolation='multilinear')

                targets.append(target)

        req = gf.Request(
            sources=[source],
            targets=targets)

        req.regularize()

        resp = self.get_engine().process(req)
        traces = resp.pyrocko_traces()

        for tr in traces:
            tr.set_ydata(num.diff(tr.ydata) / tr.deltat)

        self.add_traces(traces)

    def mechanism_from_event(self):

        event = self.get_viewer().get_active_event()

        if event is None:
            self.fail('No active event set.')

        if event.moment_tensor is not None:
            mt = event.moment_tensor.m()
        else:
            self.fail('No source mechanism available for event %s.' % event.name)
        
        self.set_parameter('lat', event.lat)
        self.set_parameter('lon', event.lon)
        self.set_parameter('depth_km', event.depth/km)

        strike, dip, slip_rake = event.moment_tensor.both_strike_dip_rake()[0]
        moment = event.moment_tensor.scalar_moment()
        self.set_parameter('magnitude', moment_tensor.moment_to_magnitude(moment))
        self.set_parameter('strike', strike)
        self.set_parameter('dip', dip)
        self.set_parameter('rake', slip_rake)

    def add_store(self):
        self._engine = self.get_engine()
        superdir = str(QFileDialog.getExistingDirectory(None,
                                 'Open working directory',
                                 '~',
                                 QFileDialog.ShowDirsOnly))
        self._engine.store_superdirs.append( superdir)
        self.store_ids = self._engine.get_store_ids()
        self.set_parameter_choices('store_id', self.store_ids)


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Seismosizer() ]


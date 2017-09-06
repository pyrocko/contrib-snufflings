import numpy as num
import math
import os

from PyQt5.QtCore import *
from pyrocko import moment_tensor, model, trace
from pyrocko.gui.snuffling import Snuffling, Param, Choice, Switch, EventMarker
from pyrocko import gf
from pyrocko.ahfullgreen import add_seismogram, Gauss, Impulse

km = 1000.


class Ahfullgreen(Snuffling):

    def __init__(self):
        Snuffling.__init__(self)

    def setup(self):
        '''Customization of the snuffling.'''
        self.set_name('Ahfullgreen')
        self.add_parameter(Param('Time', 'time', 0.0, -50., 50.))
        self.add_parameter(Param('North shift [km]', 'north_km', 10.0, -100., 100.))
        self.add_parameter(Param('East shift [km]', 'east_km', 10.0, -100., 100.))
        self.add_parameter(Param('Depth', 'depth_km', 10.0, 0.0, 600.0))
        self.add_parameter(Param('Moment', 'moment', 1., 1., 1E10))
        self.add_parameter(Param('Strike', 'strike', 0., -180., 180.))
        self.add_parameter(Param('Dip', 'dip', 90., 0., 90.))
        self.add_parameter(Param('Rake', 'rake', 0., -180., 180.))
        self.add_parameter(Param('sampling rate [Hz]', 'fsampling', 1000., 1.,
                                 10000.0))
        self.add_parameter(Param('vp [km/s]', 'vp', 6.0, 0.0, 10.0))
        self.add_parameter(Param('vs [km/s]', 'vs', 3.0, 0.0, 10.0))
        self.add_parameter(Param('Density [kg/m3]', 'density', 3000. , 0.0, 10000.0))
        self.add_parameter(Param('Qp', 'qp', 200.0, 0.0, 10000.0))
        self.add_parameter(Param('Qs', 'qs', 100.0, 0.0, 10000.0))
        self.add_parameter(Param('tau', 'tau', 0.1, 0.0, 2.0))

        self.add_parameter(Choice(
            'source shape', 'stf', 'Impulse', ['Gauss', 'Impulse']))
        self.add_parameter(Choice(
            'Waveform type', 'quantity', 'Displacement [m]',
            ['Displacement [m]', 'Velocity [m/s]', 'Acceleration [m/s2]']))
        self.add_parameter(Switch('near field', 'want_near', True))
        self.add_parameter(Switch('intermediate field', 'want_intermediate', True))
        self.add_parameter(Switch('far field', 'want_far', True))
        self.add_trigger('Set Params from Event', self.mechanism_from_event)

        self.offline_config = None

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        olat = 0.
        olon = 0.
        f = (0., 0., 0.)
        deltat = 1./self.fsampling
        if self.stf == 'Gauss':
            stf = Gauss(self.tau)
        elif self.stf == 'Impulse':
            stf = Impulse()

        viewer = self.get_viewer()
        event = viewer.get_active_event()
        if event:
            event, stations = self.get_active_event_and_stations(missing='warn')
        else:
            event = model.Event(lat=olat, lon=olon)
            stations = []

        if not stations:
            s = model.Station(lat=olat, lon=olon, station='AFG')
            stations = [s]
            viewer.add_stations(stations)

        source = gf.DCSource(
            time=event.time+self.time,
            lat=event.lat,
            lon=event.lon,
            north_shift=self.north_km*km,
            east_shift=self.east_km*km,
            depth=self.depth_km*km,
            magnitude=moment_tensor.moment_to_magnitude(self.moment),
            strike=self.strike,
            dip=self.dip,
            rake=self.rake)

        source.regularize()

        m = EventMarker(source.pyrocko_event())
        self.add_marker(m)

        targets = []

        mt = moment_tensor.MomentTensor(
            strike=source.strike,
            dip=source.dip,
            rake=source.rake,
            moment=self.moment)

        traces = []
        for station in stations:
            xyz = (self.north_km*km, self.east_km*km, self.depth_km*km)
            r = num.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
            ns = math.ceil(r/self.vs/1.6)*2
            outx = num.zeros(int(ns))
            outy = num.zeros(int(ns))
            outz = num.zeros(int(ns))
            nsl = station.nsl()
            quantity = self.quantity.split()[0].lower()
            add_seismogram(
                self.vp*km, self.vs*km, self.density, self.qp, self.qs, xyz, f,
                mt.m6(), quantity, deltat, 0., outx, outy, outz,
                stf=stf, want_near=self.want_near,
                want_intermediate=self.want_intermediate,
                want_far=self.want_far)

            for channel, out in zip('NEZ', [outx, outy, outz]):
                tr = trace.Trace('', station.station, '', channel, deltat=deltat,
                                 tmin=source.time, ydata=out)
                traces.append(tr)
        self.add_traces(traces)

    def mechanism_from_event(self):

        event = self.get_viewer().get_active_event()

        if event is None:
            self.fail('No active event set.')

        if event.moment_tensor is not None:
            strike, dip, slip_rake = event.moment_tensor.both_strike_dip_rake()[0]
            moment = event.moment_tensor.scalar_moment()
            self.set_parameter('magnitude', moment_tensor.moment_to_magnitude(moment))
            self.set_parameter('strike', strike)
            self.set_parameter('dip', dip)
            self.set_parameter('rake', slip_rake)
        else:
            self.warn('No source mechanism available for event %s. Only setting location' % event.name)

        self.set_parameter('lat', event.lat)
        self.set_parameter('lon', event.lon)
        self.set_parameter('depth_km', event.depth/km)


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ Ahfullgreen() ]


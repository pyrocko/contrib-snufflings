from pyrocko.snuffling import Snuffling, Param, Switch, Choice
from pyrocko.pile_viewer import Marker, EventMarker, PhaseMarker
from pyrocko import orthodrome

from pyrocko import cake
import os


class CakePhase(Snuffling):
    
    
    def setup(self):
        self.set_name('Cake Phase')
        
        #self._phase_names = 'PmP ~S ~P ~P(moho)s ~P(sill-top)s ~P(sill-bottom)s Pdiff'.split()
        self._phase_names = '~P Pg Sg pP P Pdiff PKP PcP PcS PKIKP pPKIKP SSP PPS SPP PSP SP PS ~PS ~SP Pn S Sn PP PPP ScS S Sdiff SS SSS PcP SKS SKIKS'.split()

        for iphase, name in enumerate(self._phase_names):
            self.add_parameter(Switch(name, 'wantphase_%i' % iphase, iphase==0)) 
        
        self.add_parameter(Param('Global shift', 'tshift', 0., -20., 20.))

        self._phases = {}
        self._model = None

    def call(self):
        '''Main work routine of the snuffling.'''
        
        self.cleanup()

        wanted = []
        for iphase, name in enumerate(self._phase_names):
            if getattr(self, 'wantphase_%i' % iphase):
                if name in self._phases:
                    phases = self._phases[name]
                else:
                    if name.startswith('~'):
                        phases = [ cake.PhaseDef(name[1:]) ]
                    else:
                        phases = cake.PhaseDef.classic(name)

                    self._phases[name] = phases
                    for pha in phases:
                        pha.name = name

                wanted.extend(phases)
        
        if not wanted:
            return

        viewer = self.get_viewer()
        pile = self.get_pile()
        event = viewer.get_active_event()

        event, stations = self.get_active_event_and_stations()
        
        if not stations:
            self.fail('No station information available.')
        
        if not self._model:
            self._model = cake.load_model()

        for station in stations:
            dist = orthodrome.distance_accurate50m(event, station)
            depth = event.depth
            if depth is None:
                depth = 0.0
            
            rays = self._model.arrivals(phases=wanted, distances=[dist*cake.m2d], zstart=depth)
            
            for ray in rays:
                time = ray.t
                name = ray.given_phase().name
                incidence_angle = ray.incidence_angle()
                takeoff_angle = ray.takeoff_angle()

                time += event.time + self.tshift
                m = PhaseMarker([ (station.network, station.station, '*', '*') ], time, time, 2, phasename=name, event=event, incidence_angle=incidence_angle, takeoff_angle=takeoff_angle)
                self.add_marker(m)

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ CakePhase() ]


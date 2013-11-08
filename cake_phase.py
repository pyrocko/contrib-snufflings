from pyrocko.snuffling import Snuffling, Param, Switch, Choice
from pyrocko.pile_viewer import PhaseMarker
from pyrocko import orthodrome

from pyrocko import cake


class CakePhase(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Theoretical Phase Arrivals</h1>
    <p>
    This snuffling uses pyrocko's <a href="http://emolch.github.io/pyrocko/v0.3/cake_doc.html">Cake</a>
    module to calculate seismic rays for layered earth models. </p>
    <p>
    <b>Parameters:</b><br />
        <b>&middot; Global shift</b>  -  Add time onset to phases. <br />
        <b>&middot; Add Model</b>  -  Add a model to drop down menu. (GUI reset required)<br />
        <b>&middot; Add Phase</b>  -  Add a phase definition. (GUI reset required)<br />
    </p>
    <p>
    Instructions and information on Cake's syntax of seismic rays can be found in the <a href="http://emolch.github.io/pyrocko/v0.3/cake_doc.html#cmdoption-cake--phase">Cake documentation</a>.
    </p>
    </body>
    </html>
    '''


    def setup(self):
        self.set_name('Cake Phase')

        #self._phase_names = 'PmP ~S ~P ~P(moho)s ~P(sill-top)s ~P(sill-bottom)s Pdiff'.split()
        self._phase_names = '~P Pg Sg pP P Pdiff PKP PcP PcS PKIKP pPKIKP SSP PPS SPP PSP SP PS ~PS ~SP Pn S Sn PP PPP ScS S Sdiff SS SSS PcP SKS SKIKS'.split()

        for iphase, name in enumerate(self._phase_names):
            self.add_parameter(Switch(name, 'wantphase_%i' % iphase, iphase==0))
        self.model_choice = Choice('Model', 'chosen_model', cake.builtin_models()[0], (cake.builtin_models()))
        self.add_parameter(self.model_choice)
        self.add_parameter(Param('Global shift', 'tshift', 0., -20., 20.))
        self.add_trigger('Add Phase', self.add_phase_definition)
        self.add_trigger('Add Model', self.add_model_to_choice)
        self._phases = {}
        self._model = None

    def call(self):

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
        if event is None:
            self.fail('No active event is marked.')

        stations = [ s for s in viewer.stations.values() if s.station in pile.stations ]

        if not stations:
            self.fail('No station information available.')

        if not self._model:
            self._model = cake.load_model(self.chosen_model)

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

    def add_model_to_choice(self):
        '''
        Called from trigger 'Add Model'.
        Adds another choice to the drop down choice menu.
        Requires a reset of the GUI.
        '''
        in_model = self.input_filename('Load Model')
        self.model_choice.choices.append(in_model)
        self.reset_gui()

    def add_phase_definition(self):
        '''
        Called from trigger 'Add Phase Definition'.
        Adds another phase option.
        Requires a reset of the GUI.
        '''
        phase_def = self.input_dialog('Add New Phase' , 'Enter Phase Definition')
        self.add_parameter(Switch(phase_def, 'wantphase_%s'%str(len(self._phases)+1), True))
        self.reset_gui()

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ CakePhase() ]


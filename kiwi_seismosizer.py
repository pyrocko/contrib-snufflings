
import re
import copy
import logging

from pyrocko.gui.snuffling import Snuffling, Param
from pyrocko import model, eventdata

logger = logging.getLogger('seismosizer')
# logger.setLevel(logging.DEBUG)

try:
    from tunguska import seismosizer, gfdb, source, glue, filtering, misfit, config
    from tunguska.phase import Taper, Timing
    config.exit_on_fatal = False
    _tunguska = True
    d2u = source.d2u
    u2d = source.u2d
    EventDataToKiwi = glue.EventDataToKiwi
except ImportError as _import_error:
    _tunguska = False
    EventDataToKiwi = object
#seismosizer.logger.setLevel(logging.DEBUG)


km = 1000.


class EventDataConverter(EventDataToKiwi):

    def __init__(self, *args, **kwargs):
        self._gfdb = kwargs.pop('gfdb')
        glue.EventDataToKiwi.__init__(self, *args, **kwargs)

    def get_preprocessed_traces(self):

        event = self._acc.get_events()[0]
        stations = self._acc.get_stations(relative_event=event)
        pile = self._acc.get_pile()

        tmin_phase = Timing('begin', -10)
        tmax_phase = Timing('end', +10)

        stations_by_distance = sorted( stations.values(), key=lambda s: s.dist_m)
        alltraces = []
        depth = event.depth
        if depth is None:
            depth = 10.*km
        for s in stations_by_distance:
            tmin, tmax = [event.time + phase(s.dist_m, depth) for phase in (tmin_phase, tmax_phase)]
            traces = pile.all(tmin=tmin, tmax=tmax, trace_selector=lambda tr: tr.nslc_id[:3] == s.nsl())
            alltraces.extend(traces)

        return alltraces

class KiwiSeismosizer(Snuffling):

    def __init__(self, sourcetype='moment_tensor'):
        Snuffling.__init__(self)
        self.sourcetype = sourcetype

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Kiwi Seismosizer (%s)' % self.sourcetype)
        if _tunguska:
            for pname in source.param_names(self.sourcetype):
                si = source.source_infos(self.sourcetype)[pname]
                self.add_parameter(Param('%s [ %s ]:' % (si.name.title(), si.unit), d2u(si.name), si.default, si.soft_min, si.soft_max))

        self.add_trigger('Configure view settings', self.compare_mode)
        self.add_trigger('Open GFDB...', self.open_gfdb)
        self.add_trigger('Set Params from Event', self.mechanism_from_event)
        self.set_live_update(False)
        self.db = None
        self.seis = None
        self.event = None

    def used_stations_with_channels(self):
        viewer = self.get_viewer()
        pile = self.get_pile()

        stations = copy.deepcopy(list(viewer.stations.values()))
        if not stations:
            logger.warn('No station information available; using dummy station at lat=2, lon=0')
            stations = [  model.Station('', 'SOUTH', '', lat=-2., lon=0.), model.Station('', 'NORTH', '', lat=2., lon=0.) ]

        stations_d = dict([ ((s.network, s.station), s) for s in stations ])

        tmin, tmax  = viewer.get_time_range()
        trs = pile.all( tmin=tmin, tmax=tmax, load_data=False, degap=False )
        if not trs:
            logger.warn('No data available; only synthetic traces will be used.')
            channel2comp = {}
            comp2channel = {}
            for cha,com in zip('ZNE', 'une'):
                for sta in stations:
                    stacha = sta.get_channel(cha)
                    if stacha is None:
                        sta.add_channel(model.Channel(cha))

                    comp2channel[sta.network, sta.station, com] = cha 

                channel2comp[cha] = com

            return stations, channel2comp, comp2channel

        stations_used = set()
        channel2comp = {}
        comp2channel = {}
        for tr in trs:
            try:
                sta = stations_d[tr.network, tr.station]
                channel2comp[tr.channel] = { 'z': 'u', 'e': 'e', 'n': 'n', 'r': 'a', 't': 'r' }[tr.channel[-1].lower()]
                stations_used.add(sta)
                comp2channel[tr.network, tr.station, channel2comp[tr.channel]] = tr.channel

                stacha = sta.get_channel(tr.channel)
                if stacha is None:
                    sta.add_channel(model.Channel(tr.channel))

            except (IndexError, KeyError):
                pass

        return stations_used, channel2comp, comp2channel

    def open_seismosizer(self):

        if self.seis:
            self.seis.close()
            self.seis = None

        stations, channel2comp, comp2channel = self.used_stations_with_channels()
        components = ''.join( sorted(list(set( k[2] for k in comp2channel.keys() ))))

        ed = eventdata.EventDataAccess(
                events=[ self.event ],
                stations=stations,
                datapile=self.get_pile() )

        if not self.db:
            self.open_gfdb()

        if not self.db:
            self.fail('No database set.')

        db = self.db

        ed_to_kiwi = EventDataConverter(ed,
                kiwi_component_map=channel2comp,
                gfdb=db,
                station_splitting=[components],
                station_filter=lambda sta: sta.dist_m < db.firstx+db.dx*(db.nx-1) )

        inner_misfit = misfit.InnerMisfitSetup(
            inner_norm = 'l1norm',
            tapers_by_set = [ Taper(timings=(Timing('begin', -10), Timing('begin', 0), Timing('end', 0), Timing('end', 10))) ],
            filters_by_set = [ filtering.Filter((0.01, 0.02, 0.5, 1))  ])

        logger.info('Starting up minimizer')
        seis = ed_to_kiwi.make_seismosizer(
                gfdb = db,
                local_interpolation = 'bilinear',
                effective_dt = db.dt*0.5)

        depth = self.event.depth
        if depth is None:
            depth = 10.*km
        inner_misfit.setup(seis, depth)

        self.comp2channel = comp2channel
        self.seis = seis

    def adjust_event(self):

        viewer = self.get_viewer()
        active_event = viewer.get_active_event()

        event = self.event

        if active_event:
            event = active_event

        if not event:
            event = model.Event(lat=0., lon=0., time=0.)
            logger.warn('No event selected; using dummy event at lat=0, lon=0.')

        if self.event is not event:
            self.event = event
            if self.seis:
                # restart seismosizer
                try:
                    self.open_seismosizer()

                except seismosizer.Fatal:
                    self.seis.close()
                    self.seis = None
                    self.fail('Minimizer crashed.')

    def get_current_source(self):

        d = {}
        for pname in source.param_names(self.sourcetype):
            si = source.source_infos(self.sourcetype)[pname]
            pname = d2u(si.name)
            d[pname] = getattr(self, pname)

        s = source.Source(self.sourcetype, **d)
        return s

    def call(self):
        '''Main work routine of the snuffling.'''
        if not _tunguska:
            self.fail('ImportError:\n%s '% _import_error)
        self.cleanup()
        self.adjust_event()

        try:
            km = 1000.
            if not self.seis:
                self.open_seismosizer()

            s = self.get_current_source()
            self.seis.set_source(s)

            snapshot = self.seis.get_receivers_snapshot(['syn'], [], 'plain')
            for rec in snapshot:
                traces = rec.get_traces()
                for tr in traces:
                    tr.set_channel(self.comp2channel[tr.network, tr.station, tr.channel])

                self.add_traces(traces)

        except seismosizer.Fatal:
            self.seis.close()
            self.seis = None
            self.fail('Minimizer crashed.')

    def open_gfdb(self):
        fn = self.input_filename()
        fn = re.sub(r'(\.\d+\.chunk|\.index)$', '', fn)

        db = gfdb.Gfdb(fn)
        if db:
            depthold = self.depth
            depthmin = db.firstz
            depthmax = db.firstz+(db.nz-1)*db.dz
            self.set_parameter_range('depth', depthmin, depthmax)
            if depthmin <= depthold <= depthmax:
                self.set_parameter('depth', depthold)
            else:
                self.set_parameter('depth', 0.5*(depthmin + depthmax))

            self.db = db

    def mechanism_from_event(self):

        self.adjust_event()

        event = self.event
        if event.moment_tensor is not None:
            mt = event.moment_tensor.m()
            logger.info('Moment tensor given with event is:\n%s' % event.moment_tensor)
        else:
            self.fail('No source mechanism available for event %s.' % event.name)

        self.set_parameter('depth', event.depth)

        if self.sourcetype in ('moment_tensor', 'mt_eikonal'):
            self.set_parameter('mxx', mt[0,0])
            self.set_parameter('myy', mt[1,1])
            self.set_parameter('mzz', mt[2,2])
            self.set_parameter('mxy', mt[0,1])
            self.set_parameter('mxz', mt[0,2])
            self.set_parameter('myz', mt[1,2])
        elif self.sourcetype in ('eikonal', 'bilateral', 'circular'):
            strike, dip, slip_rake = event.moment_tensor.both_strike_dip_rake()[0]
            moment = event.moment_tensor.scalar_moment()
            self.set_parameter('strike', strike)
            self.set_parameter('dip', dip)
            self.set_parameter('slip_rake', slip_rake)
            self.set_parameter('moment', moment)

        # set extensions to zero
        if self.sourcetype in ('eikonal', 'mt_eikonal'):
            self.get_paramater('bord_radius').set_value(0.0)
        elif self.sourcetype == 'bilateral':
            self.set_parameter('length_a', 0.0)
            self.set_parameter('length_b', 0.0)
            self.set_parameter('width', 0.0)
        elif self.sourcetype == 'circular':
            self.set_parameter('radius', 0.0)

        if event.duration is not None:
            self.set_parameter('rise_time', event.duration)
            self.set_parameter_range('rise_time', 0., 200.)

    def compare_mode(self):
        v = self.get_viewer()
        v.menuitem_colortraces.setChecked(True)
        v.menuitem_showboxes.setChecked(False)
        v.menuitems_scaling[2][0].setChecked(True)
        v.menuitems_sorting[5][0].setChecked(True)
        v.menuitems_ssorting[1][0].setChecked(True)
        v.scalingmode_change()
        v.sortingmode_change()
        v.s_sortingmode_change()
        v.update()

    def delete_gui(self):
        Snuffling.delete_gui(self)
        if self.seis:
            logger.info('Shutting down minimizer')
            self.seis.close()


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ KiwiSeismosizer('moment_tensor'), KiwiSeismosizer('bilateral') ]

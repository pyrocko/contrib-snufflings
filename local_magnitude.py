import os
import copy
import numpy as num
from collections import defaultdict
from pyrocko.snuffling import Snuffling, Param, PhaseMarker, Switch, Choice, \
    EventMarker
from pyrocko import guts, orthodrome, trace, util
from pyrocko.gui_util import to01
from pyrocko.plot import graph_colors

km = 1000.

wood_anderson_response = trace.PoleZeroResponse(
    zeros=[0., 0.],
    poles=[(-5.49779 - 5.60886j), (-5.49779 + 5.60886j)],
    constant=1.
)

class LocalMagnitudeSnuffling(Snuffling):

    '''
    Local Magnitude Estimation
    --------------------------

    Main Control High- and Lowpass filters are applied to the data before
    simulating the Wood-Anderson receiver.

    The suggested default values for the geometrical spreading, anelastic
    attenuation and static magnification are those recommended by IASPEI.

    For correct estimates these values have to be calibrated for the region
    under investigation.

    For further information:
    http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:816929:1/component/escidoc:816928/IS_3.3_rev1.pdf

    Waveform data either need to have unit 'meters' or if instument
    responses have to be removed select *Needs restitution* and use the file
    browser to read responses.

    References:
    - Bormann,  P. and Dewey J., 2012. The new IASPEI standards for determining magnitudes from digital data and their relation to classical magnitudes. doi:
    - Hutton, K.L. and Boore D.M., 1987. The ML scale in southern California. Bull. seism. Soc. Am., 77, 2074-2094
    - Richter C.F., 1935. An instrumental earthquake magnitude scale, Bull. seism. Soc. Am., 25, 1-32.
    '''
    def setup(self):
        self._responses = None
        self.add_parameter(Param(
            'geom. spreading', 'const_a', 1.11, 1., 2.))
        self.add_parameter(Param(
            'anelastic attenuation', 'const_b', 0.00189, 0., 1.))
        self.add_parameter(Param(
            'static magnification', 'const_c', -2.09, -5., 5.))
        self.add_parameter(Param(
            'Duration for "fixed" time window',
            'duration_fixed', 200., 1., 500.))
        self.add_parameter(Choice(
            'Time window', 'time_window', 'visible / selected',
            ['visible / selected', 'fixed', 'distance dependent']))
        self.add_parameter(Choice(
            'Apply to', 'apply_to', 'active event',
            ['active event', 'selected events', 'all events']))
        self.add_parameter(Switch(
            'Show restituted traces', 'show_restituded_traces', False))
        self.add_parameter(Switch(
            'Mark readings', 'show_markers', False))
        self.add_parameter(Switch(
            'Show plot', 'show_plot', False))
        self.add_parameter(Switch(
            'Show message', 'do_show_message', True))
        self.add_parameter(Switch(
            'Needs restitution', 'needs_restitution', False))
        self.add_parameter(Switch(
            'Set event magnitude', 'modify_inplace', False))

        self.set_name('Local Magnitude')
        self.set_live_update(False)

        self.vmin = 1500.
        self.vmax = 6000.

    def read_responses(self, dirname):
        responses = {}
        entries = os.listdir(dirname)
        for entry in entries:
            if entry.endswith('.pf'):
                key = tuple(entry[:-3].split('.'))
                fn = os.path.join(dirname, entry)
                resp = guts.load(filename=fn)
                responses[key] = resp

        return responses

    def local_magnitude(self, distance, amplitude):

        return num.log10(amplitude*1.0e9) + \
            self.const_a*num.log10(distance/km) + \
            self.const_b*distance/km + self.const_c

    def get_response(self, nslc):
        if self._responses is None:
            self._responses = self.read_responses(self.input_directory())

        n, s, l, c = nslc
        for k in [(n, s, l, c), (n, s, c), (s, c), (s,)]:
            if k in self._responses:
                return self._responses[k]

        self.fail(
            'no response information available for trace %s.%s.%s.%s' % nslc)

    def get_traces(self, event, stations, trace_selector, tpad):
        p = self.get_pile()
        trace_selector_viewer = self.get_viewer_trace_selector('visible')
        if self.time_window == 'distance dependant':
            for station in stations:
                distance = orthodrome.distance_accurate50m(event, station)
                tmin = distance / self.vmax
                tmax = (distance + event.depth) / self.vmin

                for trs in p.chopper(
                        tmin=event.time + tmin,
                        tmax=event.time + tmax,
                        tpad=tpad,
                        trace_selector=lambda tr: (
                            trace_selector(tr) and
                            trace_selector_viewer(tr) and
                            tr.nslc_id[:3] == station.nsl())):

                    for tr in trs:
                        yield tr

        elif self.time_window == 'fixed':
            tmin = 0.
            tmax = self.duration_fixed
            for trs in p.chopper(
                    tmin=event.time + tmin,
                    tmax=event.time + tmax,
                    tpad=tpad,
                    trace_selector=lambda tr: (
                        trace_selector(tr) and
                        trace_selector_viewer(tr))):

                for tr in trs:
                    yield tr

        else:
            for trs in self.chopper_selected_traces(
                    fallback=True, tpad=tpad,
                    trace_selector=trace_selector, mode='inview'):

                for tr in trs:
                    yield tr

    def call(self):
        self.cleanup()

        if self.apply_to == 'active event':
            event, _ = self.get_active_event_and_stations()
            events = [event]

        elif self.apply_to == 'selected events':
            events = [m.get_event() for m in self.get_selected_event_markers()]

        elif self.apply_to == 'all events':
            events = []
            for m in self.get_markers():
                if isinstance(m, EventMarker):
                    events.append(m.get_event())

        events.sort(key=lambda ev: ev.time)

        if not events:
            self.fail('no event selected')

        if self.time_window == 'visible / selected' and len(events) != 1:
            self.fail('cannot work with multiple events with "visible / '
                      'selected" time window setting')

        stations_dict = dict((s.nsl(), s) for s in self.get_stations())

        markers = []
        local_magnitudes = []
        viewer = self.get_viewer()
        fmin = viewer.highpass
        fmax = viewer.lowpass
        if not fmin or not fmax:
            self.fail('Main Controls Highpass and Lowpass filters have to be set') # noqa

        for event in events:
            mags = defaultdict(list)
            tpad = 2./fmin

            def trace_selector(tr):
                c = tr.channel.upper()
                return c.endswith('E') or c.endswith('N') or \
                    tr.location.endswith('_rest')

            distances = {}
            rest_traces = []

            event2 = copy.deepcopy(event)

            for tr in self.get_traces(
                    event, stations_dict.values(), trace_selector, tpad):

                nslc = tr.nslc_id

                try:
                    tr.highpass(4, fmin, nyquist_exception=True)
                    tr.lowpass(4, fmax, nyquist_exception=True)
                except trace.AboveNyquist, e:
                    self.fail(str(e))

                try:
                    station = stations_dict[nslc[:3]]
                except KeyError as e:
                    print e

                if self.needs_restitution:
                    resp = self.get_response(nslc)

                    try:
                        tr_vel = tr.transfer(
                            tfade=tpad,
                            freqlimits=(
                                fmin*0.5, fmin,
                                fmax, fmax*2.0),
                            transfer_function=resp,
                            invert=True)
                    except trace.TraceTooShort, e:
                        self.fail(str(e))
                        continue

                else:
                    try:
                        tr_vel = tr.transfer(
                            tfade=tpad,
                            freqlimits=(
                                fmin*0.5, fmin,
                                fmax, fmax*2.0),
                            transfer_function=wood_anderson_response,
                            invert=False)
                    except trace.TraceTooShort, e:
                        self.fail(str(e))
                        continue

                distance = orthodrome.distance_accurate50m(event, station)

                tr_vel.set_codes(location=tr_vel.location+'_rest')
                tr_vel.meta = dict(tabu=True)
                t_of_max, amplitude = tr_vel.absmax()

                if self.show_restituded_traces:
                    rest_traces.append(tr_vel)
                    m_nslc = tr_vel.nslc_id
                else:
                    m_nslc = tr.nslc_id

                mag = self.local_magnitude(distance, amplitude)
                if self.show_markers:
                    markers.append(PhaseMarker(
                        [m_nslc],
                        t_of_max, t_of_max, 1, phasename='%3.1f' % mag,
                        event=event2))

                mags[nslc[:2]].append(mag)
                distances[nslc[:2]] = distance

            if not mags:
                continue

            if rest_traces:
                self.add_traces(rest_traces)

            for k in mags:
                mags[k] = max(mags[k])

            local_magnitude = round(num.median(mags.values()), 1)

            if self.show_plot:
                data = []
                for k in mags:
                    data.append((distances[k], mags[k]))

                dists, mags_arr = num.array(data).T

                dists /= km
                fig = self.figure()
                axes = fig.add_subplot(1, 1, 1)
                axes.plot(dists, mags_arr, 'o', color=to01(graph_colors[0]))
                for x, y, label in zip(dists, mags_arr, mags.keys()):
                    axes.text(x, y, '.'.join(label))

                axes.axhline(local_magnitude, color=to01(graph_colors[0]))
                mag_std = num.std(mags.values())

                msg = 'local magnitude: %s, std: %s' % \
                    (round(local_magnitude, 1),
                        round(mag_std, 1))
                axes.text(max(dists), local_magnitude, msg,
                        verticalalignment='bottom',
                        horizontalalignment='right')

                axes.axhspan(
                        local_magnitude-mag_std,
                        local_magnitude+mag_std,
                        alpha=0.1)

                axes.set_xlabel('Distance [km]')
                axes.set_ylabel('Local Magnitude')
                fig.canvas.draw()

            local_magnitudes.append(local_magnitude)

        if self.modify_inplace:
            event.magnitude = local_magnitude

        if markers:
            self.add_markers(markers)

        if not local_magnitudes:
            self.fail('no results')

        if self.do_show_message:
            self.show_message('Magnitude', ', '.join(
                '%3.1f' % mag for mag in local_magnitudes))


def __snufflings__():
    return [LocalMagnitudeSnuffling()]

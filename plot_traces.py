import numpy as num
from pyrocko.gui.snuffling import Snuffling, Param, Choice, Switch
from pyrocko import orthodrome as ortho
from pyrocko import util


def station_key(tr):
    return (tr.network, tr.station, tr.location)


class TracePlotter(Snuffling):
    '''
    <html>
    <body>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Plot Traces with Reduced Velocity</h1>
    <p>
    Use the <b>Reduction Velocity</b> to shift traces dependent on epicentral
    distance to the activated event.
    If <b>Auto-Run</b> is activated the figure is updated automatically when
    modifying a value on the panel.<br>
    When <b>saving a figure</b> accepted file endings are 'eps', 'png', 'jpeg',
    'pdf', and 'svg'.
    </body>
    </html>
    '''
    def setup(self):
        self.set_name("Plot Waveforms")
        self.add_parameter(
            Switch('Include Selected Markers', 'add_markers',  False))

        self.add_parameter(Switch('Fill positive', 'fill_between',  False))
        self.add_parameter(
            Param(
                'Reduction Velocity [km/s]',
                't_red', 20., 1., 20., high_is_none=True))

        self.add_parameter(Param(
            'Amplitude Gain', 'yscale', 1., 0.1, 100.))
        self.add_parameter(Choice(
            'Pre-scale Amplitudes', 'ampl_scaler', 'trace min/max',
            ['total min/max', 'trace min/max', 'standard deviation']))
        self.add_trigger('Save Last Figure', self.save)
        self.set_live_update(False)
        self.fig = None

    def call(self):

        self.cleanup()
        viewer = self.get_viewer()

        vtmin, vtmax = viewer.get_time_range()
        pile = self.get_pile()
        traces = [
            tr for tr in pile.chopper(
                tmin=vtmin, tmax=vtmax, trace_selector=viewer.trace_selector)]

        event, stations = self.get_active_event_and_stations()
        traces = [tr for trs in traces for tr in trs]

        trace_nsls = {station_key(tr) for tr in traces}
        stations = [s for s in self.get_stations() if
                    s.nsl() in trace_nsls]

        distances = [
            ortho.distance_accurate50m(event, s)/1000. for s in stations]
        distances = dict(zip([s.nsl() for s in stations], distances))

        matching_traces = [x for x in traces if util.match_nslc(
                            self.get_station_patterns(stations), x.nslc_id)]
        if self.add_markers:
            markers = self.get_markers()
            markers = [
                m for m in markers if m.tmax <= vtmax and
                m.tmin >= vtmin and m.selected]

            markers = dict(zip([tuple(m.nslc_ids) for m in markers], markers))

        if self.fig is None or self.fframe.closed or not self._live_update:
            self.fframe = self.pylab(get='figure_frame')
            self.fig = self.fframe.gcf()

        if self._live_update:
            self.fig.clf()

        maxd = max(distances.values())
        mind = min(distances.values())

        ymin = mind-0.06*(maxd-mind)
        ymax = maxd+0.06*(maxd-mind)
        ax = self.fig.add_subplot(111)
        xmin = 9E9
        xmax = -xmin
        texts = []
        manual_scale = 0.1 * (maxd-mind)*self.yscale

        if self.ampl_scaler == 'total min/max':
            max_trace = max(
                matching_traces, key=lambda x: max(abs(x.get_ydata())))

            tr_maxy = max(abs(max_trace.get_ydata()))
            ampl_scale = float(tr_maxy)

        for tr in matching_traces:
            if viewer.highpass:
                tr.highpass(4, viewer.highpass)
            if viewer.lowpass:
                tr.lowpass(4, viewer.lowpass)
            if tr.nslc_id[:3] not in distances.keys():
                continue

            if self.t_red:
                red = distances[tr.nslc_id[:3]]/self.t_red
            else:
                red = 0.
            y_pos = distances[tr.nslc_id[:3]]
            xdata = tr.get_xdata()-red-event.time
            xmin = min(xmin, min(xdata))
            xmax = max(xmax, max(xdata))
            tr_ydata = tr.get_ydata()
            if self.ampl_scaler == 'trace min/max':
                ampl_scale = float(max(abs(tr_ydata)))
            elif self.ampl_scaler == 'standard deviation':
                ampl_scale = float(num.std(tr_ydata))
            ydata = (tr_ydata/ampl_scale * manual_scale) + y_pos
            ax.plot(xdata, ydata, c='black', linewidth=0.2)

            if self.fill_between:
                ax.fill_between(
                    xdata, y_pos, ydata, where=ydata > y_pos, color='black',
                    alpha=0.5)

            texts.append(
                ax.text(
                    xmax, y_pos, '%s.%s.%s.%s' % tr.nslc_id,
                    horizontalalignment='right', fontsize=6.))

            if self.add_markers:
                for ids, m in markers.items():
                    if m.match_nslc(tr.nslc_id) or ids == ():
                        c = m.select_color(m.color_b)
                        c = [ci/255. for ci in c]
                        t = m.tmin
                        x = [t-red-event.time, t-red-event.time]
                        y = [y_pos-(maxd-mind)*0.025, y_pos+(maxd-mind)*0.025]
                        ax.plot(x, y, linewidth=1, color=c)
                        label = m.get_label()
                        if not label:
                            label = ''

                        ax.text(x[1]-x[1]*0.005, y[1], label, color=c,
                                fontsize=6,
                                verticalalignment='top',
                                horizontalalignment='right')

        for txt in texts:
            txt.set_x(xmax)

        vred_str = '= '+str(round(self.t_red, 2)) + 'km/s' if self.t_red \
            else 'off'

        ax.text(0.5, 0.01, 'time window: %s - %s  |   Reduction velocity %s' %
                (util.tts(vtmin), util.tts(vtmax), vred_str),
                verticalalignment='bottom', horizontalalignment='center',
                transform=self.fig.transFigure)

        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])
        ax.set_ylabel('Distance [km]')
        ax.set_xlabel('(red.) Time [s]')
        self.fig.canvas.draw()

    def save(self):
        fn = self.output_filename('Select Filename', 'snuffled_traces.png')
        self.fig.savefig(fn, pad_inches=0.1, bbox_inches='tight', dpi=320)

    def set_center_latlon(self):
        self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(
                self.get_stations())

        self.set_parameter('lat_c', self.lat_c)
        self.set_parameter('lon_c', self.lon_c)

    def get_station_patterns(self, stations):
        return ['%s.%s.%s.*' % s.nsl() for s in stations]


def __snufflings__():
    return [TracePlotter()]

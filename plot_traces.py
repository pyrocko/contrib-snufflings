from pyrocko.snuffling import Snuffling, Param
from pyrocko import orthodrome as ortho
from pyrocko import util


class TracePlotter(Snuffling):
    '''
    <html>
    <body>
    <h1>Trace Plotter</h1>
    </body>
    </html>
    '''
    def setup(self):
        self.set_name("Plot Waveforms")
        self.add_parameter(Param('Reduction velocity[km/s]', 't_red', 1., 1., 10., low_is_none=True))
        self.set_live_update(False)

    def call(self):

        self.cleanup()

        event, stations = self.get_active_event_and_stations()
        distances = [ortho.distance_accurate50m(event, s) for s in stations]

        maxd = max(distances)
        mind = min(distances)
        distances = dict(zip([s.nsl() for s in stations], distances))
        channels = set()
        traces = list(self.chopper_selected_traces(fallback=True))
        traces = [tr for trs in traces for tr in trs ]
        matching_traces = filter(lambda x: util.match_nslc(
                            self.get_station_patterns(stations), x.nslc_id), traces)
        fig = self.pylab(get='figure')
        ax = None
        viewer = self.get_viewer()
        axs = {}
        for tr in matching_traces:
            if viewer.highpass:
                tr.highpass(4, viewer.highpass)
            if viewer.lowpass:
                tr.lowpass(4, viewer.lowpass)
            if tr.nslc_id[:3] not in distances.keys():
                continue

            ax = fig.add_axes([0.1, 0.1+0.9*((distances[tr.nslc_id[:3]]-mind)/(maxd-mind)), 0.93, 1./len(matching_traces)],
                              sharex=ax)
            if self.t_red:
                red = distances[tr.nslc_id[:3]]/self.t_red/1000.
            else:
                red = 0.
            ax.plot(tr.get_xdata()-red-event.time,
                    tr.get_ydata(),
                    c='black',
                    linewidth=0.1)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.patch.set_visible(False)
            ax.text(
                0., 0.5, '%s.%s.%s.%s' % tr.nslc_id, transform=ax.transAxes, 
                horizontalalignment='right', fontsize=4.)
            for item in ax.spines.values():
                item.set_visible(False)
            axs[distances[tr.nslc_id[:3]]] = ax
        axs[mind].get_xaxis().set_visible(True)
        fig.canvas.draw()

    def save(self):
        pass

    def set_center_latlon(self):
        self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(self.get_stations())
        self.set_parameter('lat_c', self.lat_c)
        self.set_parameter('lon_c', self.lon_c)

    def get_station_patterns(self, stations):
        return ['%s.%s.%s.*' % s.nsl() for s in stations]

def __snufflings__():
    return [TracePlotter()]

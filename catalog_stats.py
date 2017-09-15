from pyrocko.gui.snuffling import Snuffling, Param, Choice
from pyrocko.gui.util import EventMarker
from pyrocko import util
from pyrocko import moment_tensor
from pyrocko.orthodrome import distance_accurate50m as distance
import matplotlib.gridspec as gridspec
import matplotlib.dates as md
import datetime
import numpy as num


class CumEvent(Snuffling):

    '''
    <html>
    <body>
    <h1>Cumulative Number of Events</h1>
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Catalog Statistics')
        self.add_parameter(
            Param('Latitude:', 'lat', 90, -90., 90., high_is_none=True))
        self.add_parameter(
            Param('Longitude:', 'lon', 180., -180, 180., high_is_none=True))
        self.add_parameter(
            Param('Maximum Distance [km]:', 'maxd', 20000., 0., 20000.,
                  high_is_none=True))
        self.add_parameter(Choice('Event distribution', 'variation',
                                  'daily', ['daily', 'annual']))
        self.add_trigger('Save Figure', self.save_as)
        self.set_live_update(False)
        self.fig = None
        self.cli_mode = False

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        viewer = self.get_viewer()
        tmin, tmax = self.get_selected_time_range(fallback=True)
        event_markers = [x for x in viewer.markers if x.tmin >= tmin and
                         x.tmax <= tmax]

        event_markers = [x for x in event_markers if isinstance(x, EventMarker)]
        if self.maxd:
            event_markers = [x for x in event_markers if
                             distance(self, x._event) <= self.maxd*1000.]

        if event_markers == []:
            self.fail('No events in selected area found')
        events = [m.get_event() for m in event_markers]

        self.make_time_line(events)

    def make_time_line(self, events):

        if self.cli_mode:
            import matplotlib.pyplot as plt
            self.fig = plt.figure()
        else:
            fframe = self.figure_frame()
            self.fig = fframe.gcf()

        gs = gridspec.GridSpec(2, 1)
        gs.update(hspace=0.005, top=0.95)
        gs1 = gridspec.GridSpec(2, 1)
        gs1.update(bottom=0.06, hspace=0.01)

        ax = self.fig.add_subplot(gs[0])
        ax1 = ax.twinx()
        ax2 = self.fig.add_subplot(gs1[-1])
        events.sort(key=lambda x: x.time)
        magnitudes = []
        cum_events = num.cumsum(num.ones(len(events)))

        for e in events:
            if e.moment_tensor is not None:
                magnitudes.append(e.moment_tensor.magnitude)
            else:
                magnitudes.append(e.magnitude)
            if magnitudes[-1] is None:
                magnitudes.pop()
                magnitudes.append(0.)

        magnitudes = moment_tensor.magnitude_to_moment(num.array(magnitudes))
        cum_events_magnitude = num.cumsum(magnitudes)
        times = num.array([e.time for e in events])
        timeslabels = [datetime.datetime(1970, 1, 1) +
                       datetime.timedelta(seconds=t) for t in times]

        ax.plot(timeslabels, cum_events)
        ax.set_ylabel('Cumulative number of events')
        ax.axes.get_xaxis().set_ticklabels([])
        ax1.plot(timeslabels, cum_events_magnitude, '-b')
        ax1.set_ylabel('Cumulative moment [Nm]')
        xfmt = md.DateFormatter('%Y-%m-%d')
        ax1.xaxis.set_major_formatter(xfmt)
        ax.yaxis.label.set_color('r')
        ax1.yaxis.label.set_color('b')

        ax.grid(True)
        ax1.grid(True)

        t0 = min(times)
        if self.variation == 'daily':
            normalization = 60 * 60
            nbins = 24
            t0_day = util.day_start(t0)
            ax2.set_xlabel('hour of day')

        elif self.variation == 'annual':
            normalization = 60 * 60 * 24
            nbins = 365
            t0_day = util.year_start(t0)
            ax2.set_xlabel('day of year')

        binned = (times - t0_day) % (nbins * normalization)
        ax2.hist(binned/normalization, bins=nbins, color='grey')
        ax2.set_ylabel('Number of events')
        ax2.set_xlim((0, nbins))

        if self.cli_mode:
            plt.show()
        else:
            self.fig.canvas.draw()

    def save_as(self):
        if self.fig:
            fn = self.output_filename()
            self.fig.savefig(fn,
                             pad_inches=0.05,
                             bbox_inches='tight')

    def configure_cli_parser(self, parser):
        parser.add_option(
            '--events',
            dest='events_filename',
            default=None,
            metavar='FILENAME',
            help='Read events from FILENAME')


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''

    return [CumEvent()]

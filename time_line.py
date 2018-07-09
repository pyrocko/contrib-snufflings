from pyrocko.gui.snuffling import Snuffling, Param, Choice, Switch
from pyrocko.gui.util import EventMarker
from pyrocko import orthodrome
from pyrocko import moment_tensor
from pyrocko.orthodrome import distance_accurate50m as distance
from pyrocko import util, model
import matplotlib.dates as mdates
from matplotlib import cm
import numpy as num
from datetime import datetime


km = 1000.
cmap = cm.RdYlBu
cmaps = {'Red-Yellow-Blue': 'RdYlBu',
         'Viridis': 'viridis',
         'Magma': 'magma'}
save_cmaps = {}
for k, v in cmaps.items():
    try:
         x = getattr(cm, v)
         save_cmaps[k] = x
    except AttributeError:
         continue


class TimeLine(Snuffling):

    '''
    <html>
    <body>
    <h1>Temporal Evolution of Seismicity</h1>

    The considered region can be limited by defining one central coordinate
    and a maximum epicentral distance.
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Time Line')
        self.add_parameter(
            Param('Latitude:', 'lat', 90, -90., 90., high_is_none=True))
        self.add_parameter(
            Param('Longitude:', 'lon', 180., -180, 180., high_is_none=True))
        self.add_parameter(
            Param('Maximum Distance [km]:', 'maxd', 20000., 0., 20000.,
                  high_is_none=True))
        self.add_parameter(
            Choice('Color by', 'color_by', 'time',
                   ['time', 'longitude', 'latitude', 'magnitude', 'depth', 'kind']))
        self.add_parameter(Choice('Colormap', 'cmap_selector',
                                  'Red-Yellow-Blue', list(save_cmaps.keys())))
        self.add_parameter(Choice('Coordinate system', 'coord_system',
                                  'Lat/Lon', ['Lat/Lon', 'cartesian']))
        self.add_parameter(Switch('Show stations', 'show_stations', False))
        self.add_trigger('Save Figure', self.save_as)
        self.set_live_update(False)
        self.fig = None
        self.cli_mode = False

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        cmap = save_cmaps[self.cmap_selector]
        viewer = self.get_viewer()
        tmin, tmax = self.get_selected_time_range(fallback=True)
        event_markers = filter(lambda x: x.tmin >= tmin and x.tmax <= tmax,
                               viewer.markers)

        event_markers = filter(
            lambda x: isinstance(x, EventMarker), event_markers)

        if self.maxd:
            event_markers = filter(
                lambda x: distance(self, x._event) <= self.maxd*km,
                event_markers)

        if event_markers == []:
            self.fail('No events in selected area found')

        event_markers = list(event_markers)
        stations = self.get_stations() if self.show_stations else None
        self.make_time_line(event_markers, stations, cmap=cmap)

    def make_time_line(self, markers, stations=None, cmap=cmap):
        events = [m.get_event() for m in markers]
        kinds = num.array([m.kind for m in markers])
        if self.cli_mode:
            self.fig = plt.figure()
        else:
            fframe = self.figure_frame()
            self.fig = fframe.gcf()
        ax = self.fig.add_subplot(311)
        ax_cum = ax.twinx()
        ax1 = self.fig.add_subplot(323)
        ax2 = self.fig.add_subplot(325, sharex=ax1)
        ax3 = self.fig.add_subplot(324, sharey=ax1)

        num_events = len(events)
        data = num.zeros((num_events, 6))
        column_to_index = dict(zip(['magnitude', 'latitude', 'longitude', 'depth', 'time', 'kind'],
                           range(6)))
        c2i = column_to_index
        for i, e in enumerate(events):
            if e.magnitude:
                mag = e.magnitude
            else:
                mag = 0.
            data[i, :] = mag, e.lat, e.lon, e.depth, e.time, kinds[i]

        s_coords = num.array([])
        s_labels = []
        if stations is not None:
            s_coords = num.array([(s.lon, s.lat, s.elevation-s.depth) for s in stations])
            s_labels = ['.'.join(s.nsl()) for s in stations]

        isorted = num.argsort(data[:, c2i['time']])
        data = data[isorted]

        def _D(key):
            return data[:, c2i[key]]

        tmin = _D('time').min()
        tmax = _D('time').max()
        lon_max = _D('longitude').max()
        lon_min = _D('longitude').min()
        lat_max = _D('latitude').max()
        lat_min = _D('latitude').min()
        depths_min = _D('depth').min()
        depths_max = _D('depth').max()
        mags_min = _D('magnitude').min()
        mags_max = _D('magnitude').max()
        moments = moment_tensor.magnitude_to_moment(_D('magnitude'))
        dates = list(map(datetime.fromtimestamp, _D('time')))

        fds = mdates.date2num(dates)
        tday = 3600*24
        tweek = tday*7
        if tmax-tmin < 1*tday:
            hfmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        elif tmax-tmin < tweek*52:
            hfmt = mdates.DateFormatter('%Y-%m-%d')
        else:
            hfmt = mdates.DateFormatter('%Y/%m')

        color_values = _D(self.color_by)
        color_args = dict(c=color_values, vmin=color_values.min(),
                    vmax=color_values.max(), cmap=cmap)

        ax.scatter(fds, _D('magnitude'), s=20, **color_args)

        ax.xaxis.set_major_formatter(hfmt)
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.set_ylim((mags_min, mags_max*1.10))
        ax.set_xlim(map(datetime.fromtimestamp, (tmin, tmax)))
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel('Magnitude')
        init_pos = ax.get_position()

        ax_cum.plot(fds, num.cumsum(moments), 'grey')
        ax_cum.xaxis.set_major_formatter(hfmt)
        ax_cum.spines['top'].set_color('none')
        ax_cum.spines['right'].set_color('grey')
        ax_cum.set_ylabel('Cumulative seismic moment')

        lats_min = num.array([lat_min for x in range(num_events)])
        lons_min = num.array([lon_min for x in range(num_events)])

        if self.coord_system == 'cartesian':
            lats, lons = orthodrome.latlon_to_ne_numpy(
                lats_min, lons_min, _D('latitude'), _D('longitude'))

            _x = num.empty((len(s_coords), 3))
            for i, (slon, slat, sele) in enumerate(s_coords):
                n, e = orthodrome.latlon_to_ne(lat_min, lon_min, slat, slon)
                _x[i, :] = (e, n, sele)
            s_coords = _x
        else:
            lats = _D('latitude')
            lons = _D('longitude')

        s_coords = s_coords.T

        ax1.scatter(lons, lats, s=20, **color_args)
        ax1.set_aspect('equal')
        ax1.grid(True, which='both')
        ax1.set_ylabel('Northing [m]')
        ax1.get_yaxis().tick_left()

        if len(s_coords):
            ax1.scatter(s_coords[0], s_coords[1], marker='v', s=40, color='black')
            for c, sl in zip(s_coords.T, s_labels):
                ax1.text(c[0], c[1], sl, color='black')

        # bottom left plot
        ax2.scatter(lons, _D('depth'), s=20, **color_args)
        ax2.grid(True)
        ax2.set_xlabel('Easting [m]')
        ax2.set_ylabel('Depth [m]')
        ax2.get_yaxis().tick_left()
        ax2.get_xaxis().tick_bottom()
        ax2.invert_yaxis()

        ax2.text(1.1, 0, 'Origin at:\nlat=%1.3f, lon=%1.3f' %
                 (lat_min, lon_min), transform=ax2.transAxes)

        # top right plot
        ax3.scatter(_D('depth'), lats, s=20, **color_args)
        ax3.set_xlim((depths_min, depths_max))
        ax3.grid(True)
        ax3.set_xlabel('Depth [m]')
        ax3.get_xaxis().tick_bottom()
        ax3.get_yaxis().tick_right()

        self.fig.subplots_adjust(
            bottom=0.05, right=0.95, left=0.075, top=0.95, wspace=0.02, hspace=0.02)
        init_pos.y0 += 0.05
        ax.set_position(init_pos)
        ax_cum.set_position(init_pos)
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

    return [TimeLine()]

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    util.setup_logging('time_line.py', 'info')
    s = TimeLine()
    options, args, parser = s.setup_cli()

    s.cli_mode = True
    if options.events_filename:
        s.make_time_line(
            list(model.Event.load_catalog(options.events_filename)))

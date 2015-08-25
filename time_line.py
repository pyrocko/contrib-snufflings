from pyrocko.snuffling import Snuffling, Param, Choice, Switch
from pyrocko.gui_util import EventMarker
from pyrocko.orthodrome import distance_accurate50m as distance
from pyrocko import util, model
import matplotlib.dates as mdates
from matplotlib import cm
import numpy as num
from datetime import datetime

cmap = cm.jet
km = 1000.

class TimeLine(Snuffling):
    
    '''
    <html>
    <body>
    <h1>Temporal Evolution of Seismicity</h1>

    The considered region can be limited by defining one central coordinate and a
    maximum epicentral distance.
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Time Line')
        self.add_parameter(Param('Latitude:', 'lat', 90, -90., 90., high_is_none=True))
        self.add_parameter(Param('Longitude:', 'lon', 180., -180, 180., high_is_none=True))
        self.add_parameter(Param('Maximum Distance [km]:', 'maxd', 20000., 0., 20000., high_is_none=True))

        self.add_trigger('Save Figure', self.save_as)
        self.set_live_update(False)
        self.fig = None
        self.cli_mode = False

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        viewer = self.get_viewer()  
        tmin, tmax = self.get_selected_time_range(fallback=True)
        event_markers = filter(lambda x: x.tmin>=tmin and x.tmax<=tmax,
                               viewer.markers)

        event_markers = filter(lambda x: isinstance(x, EventMarker), event_markers)
        if self.maxd:
            event_markers = filter(lambda x: distance(self, x._event)<=self.maxd*km,
                               event_markers)

        if event_markers==[]:
            self.fail('No events in selected area found')
        events = [m.get_event() for m in event_markers]

        self.make_time_line(events)

    def make_time_line(self, events):

        if self.cli_mode:
            self.fig = plt.figure()
        else:
            fframe = self.figure_frame()
            self.fig = fframe.gcf()
        ax = self.fig.add_subplot(311)
        ax1 = self.fig.add_subplot(323)
        ax2 = self.fig.add_subplot(325)
        ax3 = self.fig.add_subplot(324)
        
        num_events = len(events)
        magnitudes = num.zeros(num_events)
        times = num.zeros(num_events)
        lats = num.zeros(num_events)
        lons = num.zeros(num_events)
        depths = num.zeros(num_events)
        for i, e in enumerate(events):
            if e.magnitude:
                mag = e.magnitude 
            else:
                mag = 0.
            magnitudes[i] = mag
            lats[i] = e.lat
            lons[i] = e.lon
            depths[i] = e.depth
            times[i] = e.time
        
        tmin = min(times)
        tmax = max(times)
        lon_max = lons.max()
        lon_min = lons.min()
        lat_max = lats.max()
        lat_min = lats.min()
        depths_min = depths.min()
        depths_max = depths.max()
        mags_min = magnitudes.min()
        mags_max = magnitudes.max()
        dates = map(datetime.fromtimestamp, times)

        fds = mdates.date2num(dates)
        tday = 3600*24
        tweek = tday*7
        if tmax-tmin<1*tday:
            hfmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        elif tmax-tmin<tweek*52:
            hfmt = mdates.DateFormatter('%Y-%m-%d')
        else:
            hfmt = mdates.DateFormatter('%Y/%m')

        ax.scatter(fds, 
                   magnitudes, 
                   s=20, 
                   c=times, 
                   vmin=tmin, 
                   vmax=tmax, 
                   cmap=cmap)

        ax.xaxis.set_major_formatter(hfmt)
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.set_ylim((mags_min, mags_max*1.10))
        ax.set_xlim(map(datetime.fromtimestamp, (tmin, tmax)))
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel('Magnitude')
        init_pos = ax.get_position()

        # top left plot
        ax1.scatter(lons, lats, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax1.set_xlim((lon_min, lon_max))
        ax1.set_ylim((lat_min, lat_max))
        ax1.grid(True, which='both')
        ax1.set_xticklabels([])
        ax1.set_ylabel('Lat')
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()

        # bottom left plot
        ax2.scatter(lons, depths, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax2.set_xlim((lon_min, lon_max))
        ax2.set_ylim((depths_min, depths_max))
        ax2.grid(True)
        ax2.set_xlabel('Lon')
        ax2.set_ylabel('Depth')
        ax2.get_yaxis().tick_left()
        ax2.invert_yaxis()

        # top left plot
        ax3.scatter(depths, lats, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax3.set_xlim((depths_min, depths_max))
        ax3.grid(True)
        ax3.set_ylim((lat_min, lat_max))
        ax3.set_xlabel('Depth')
        ax3.get_xaxis().tick_bottom()
        ax3.get_yaxis().tick_right()

        self.fig.subplots_adjust(bottom=0.1, 
                            right=0.9, 
                            top=0.95,
                            wspace=0.02,
                            hspace=0.02)
        init_pos.y0+=0.05
        ax.set_position(init_pos)
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
    
    return [ TimeLine() ]

if __name__=='__main__':
    import matplotlib.pyplot as plt
    util.setup_logging('time_line.py', 'info')
    s = TimeLine()
    options, args, parser = s.setup_cli()

    s.cli_mode = True
    if options.events_filename:
        s.make_time_line(list(model.Event.load_catalog(options.events_filename)))


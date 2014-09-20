from pyrocko.snuffling import Snuffling, Param, Choice, Switch
from pyrocko.gui_util import EventMarker
from pyrocko.orthodrome import distance_accurate50m as distance
from pyrocko.util import str_to_time, time_to_str
import matplotlib.dates as mdates
from matplotlib import cm
import numpy as num

cmap = cm.jet
km = 1000.

class TimeLine(Snuffling):
    
    '''
    <html>
    <body>
    <h1>Temporal Evolution of Seismicity</h1>
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Time Line')
        self.add_parameter(Param('Latitude:', 'lat', 90, -90., 90., high_is_none=True))
        self.add_parameter(Param('Longitude:', 'lon', 180., -180, 180., high_is_none=True))
        self.add_parameter(Param('Maximum Distance [km]:', 'maxd', 20000., 0., 20000., high_is_none=True))

        self.add_parameter(Switch('Save figure', 'save', False))
        self.set_live_update(False)
        
    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        viewer = self.get_viewer()
        tmin, tmax = self.get_selected_time_range(fallback=True)
        event_markers = filter(lambda x: x.tmin>tmin and x.tmax<tmax,
                               viewer.markers)
        event_markers = filter(lambda x: isinstance(x, EventMarker), event_markers)
        if self.maxd:
            event_markers = filter(lambda x: distance(self, x._event)<=self.maxd*km,
                               event_markers)

        fframe = self.figure_frame()
        fig = fframe.gcf()
        ax = fig.add_subplot(311)
        ax1 = fig.add_subplot(323)
        ax2 = fig.add_subplot(325)
        ax3 = fig.add_subplot(324)
        ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        
        num_events = len(event_markers)
        magnitudes = num.zeros(num_events)
        times = num.zeros(num_events)
        lats = num.zeros(num_events)
        lons = num.zeros(num_events)
        depths = num.zeros(num_events)
        for i, m in enumerate(event_markers):
            e = m._event
            if e.magnitude:
                mag = e.magnitude 
            else:
                mag = 0.
            magnitudes[i] = mag
            lats[i] = e.lat
            lons[i] = e.lon
            depths[i] = e.depth
            times[i] = e.time
            time_str = time_to_str(m.tmin)
            time_str = time_str.replace(',','.')
        
        lon_max = lons.max()
        lon_min = lons.min()
        lat_max = lats.max()
        lat_min = lats.min()
        depths_min = depths.min()
        depths_max = depths.max()
        mags_min = magnitudes.min()
        mags_max = magnitudes.max()

        ax.scatter(times, 
                   magnitudes, 
                   s=20, 
                   c=times, 
                   vmin=tmin, 
                   vmax=tmax, 
                   cmap=cmap)

        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.set_ylim((mags_min, mags_max*1.10))
        ax.set_xlim((tmin, tmax))
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel('Magnitude')
        init_pos = ax.get_position()

        #fig.autofmt_xdate()
        # top left plot
        ax1.scatter(lons, lats, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax1.set_xlim((lon_min, lon_max))
        ax1.set_ylim((lat_min, lat_max))
        ax1.grid(True, which='both')
        ax1.set_xticklabels([])
        ax1.set_ylabel('latitude')

        # bottom left plot
        ax2.scatter(lons, depths, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax2.set_xlim((lon_min, lon_max))
        ax2.set_ylim((depths_min, depths_max))
        ax2.grid(True)
        ax2.set_xlabel('longitude')
        ax2.set_ylabel('Depth')

        ax3.scatter(depths, lats, s=20, c=times, vmin=tmin, vmax=tmax, cmap=cmap)
        ax3.set_xlim((depths_min, depths_max))
        ax3.grid(True)
        ax3.set_ylim((lat_min, lat_max))
        ax3.yaxis.tick_right()

        fig.subplots_adjust(bottom=0.1, 
                            right=0.8, 
                            top=0.9,
                            wspace=0.0,
                            hspace=0.0)

        ax.set_position(init_pos)
        fig.canvas.draw()

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ TimeLine() ]


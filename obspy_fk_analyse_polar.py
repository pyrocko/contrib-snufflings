from pyrocko.snuffling import Param, Snuffling, pile , Choice
from pyrocko import trace
from obspy.core import UTCDateTime, stream
from obspy.core import trace as obspy_trace
from obspy.signal import array_analysis

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

def p2o_trace(ptrace, station):
    '''Convert Pyrocko trace to ObsPy trace.'''

    otr = obspy_trace.Trace(
            data = ptrace.get_ydata(),
            header=dict(
                network = ptrace.network,
                station = ptrace.station,
                location = ptrace.location,
                channel = ptrace.channel,
                delta = ptrace.deltat,
                starttime = UTCDateTime(ptrace.tmin),
                coordinates = dict(
                    latitude = station.lat,
                    longitude = station.lon,
                    elevation = station.elevation/1000. )))

    return otr

class FK(Snuffling):
    '''
    <html>
    <body>
    <h1>FK ANALYSIS</h1>
    <h2>based in ObsPy</h2>
    <p>
    <u>Prerequisites</u><br>
    This snuffling requires the ObsPy package which can be found at 
        <a href='https://github.com/obspy/obspy/wiki'> ObsPy's github page </a>
    </p>
    <p>
    The easiest way to install this package is to do:
        <PRE>
        git clone git://github.com/obspy/obspy.git
        cd obspy 
        sudo python setup.py install 
        </PRE>
    <p>
    <u>Usage</u><br>
    
     - Snuffler needs to be invoked with station information, as for example:
            snuffler testdata.mseed --stations=testdata.stations<br>

     - Zoom into the data until you see only data you desire to analyse.<br>

     - Press the 'Run' button. <br>

    The slowness is given in s/km.
    The location of the geometrical center is printed to the terminal.
    </p>
    <p>
    <u>Error Handling</u><br>
    
    If you get an error like
        <PRE>ValueError: shape mismatch: objects cannot be broadcast to a single shape</PRE>
    you need to reduce the 'Length of Sliding Window [s]' by moving the slider to the left.
    </p>
    </body>
    </html>

    '''
    
    def setup(self):
        self.set_name('FK Analysis')
        self.add_parameter(Param('Slowness range[+-]','smax',0.5,0.1, 3))
        self.add_parameter(Param('Number of slowness divisions','divisor',20,10,50))
        self.add_parameter(Param('Number of radial sections','numberOfFraction',32,4,50))
        self.add_parameter(Param('Length of Sliding Window [s]','window_lenth',1.,0.5,5.))
        #self.add_parameter(Choice('Units: ','unit','[s/km]',('[s/km]','[s/deg]')))
        self.set_live_update(False)

    def call(self):
        self.cleanup()
        viewer = self.get_viewer()

        if viewer.lowpass is None or viewer.highpass is None:
            self.fail('highpass and lowpass in viewer must be set!')

        traces = []
        for trs in self.chopper_selected_traces(fallback=True):
            for tr in trs:
                tr.downsample_to(1/20.)
                tr.lowpass(4, viewer.lowpass)
                tr.highpass(4, viewer.highpass)

            traces.extend(trs)

        if not traces:
            self.fail('no traces selected')

        tmin, tmax = trace.minmaxtime(traces, key=lambda x: None)[None]

        try:
            obspy_traces = [ p2o_trace(tr, viewer.get_station(viewer.station_key(tr)) ) for tr in traces ]

        except KeyError:
            self.fail('station information missing')

        st = stream.Stream(traces=obspy_traces)
        center = array_analysis.get_geometry(st, return_center=True)
        center_lon, center_lat, center_ele = center[len(center)-1]


        # Execute sonic
        kwargs = dict(
            # slowness grid: X min, X max, Y min, Y max, Slow Step
            sll_x=-self.smax, slm_x=self.smax, sll_y=-self.smax, slm_y=self.smax, sl_s=self.smax/self.divisor,
            # sliding window properties
            win_len=self.window_lenth, win_frac=0.1,
            # frequency properties
            frqlow=viewer.highpass, frqhigh=viewer.lowpass, prewhiten=0,
            # restrict output
            semb_thres=-1.0e9, vel_thres=-1.0e9, verbose=True, timestamp='mlabday',
            stime=UTCDateTime(tmin), etime=UTCDateTime(tmax)
        )
        

        try:
            out = array_analysis.array_processing(st, **kwargs)
        except AttributeError:
            from obspy.signal.array_analysis import sonic
            out = sonic(st, **kwargs)

        cmap = cm.hot_r
        pi = np.pi

        # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = out.T
        
        baz[baz < 0.0] += 360.

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = int(self.numberOfFraction)
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0., self.smax, N + 1)

        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = np.histogram2d(baz, slow,
                bins=[abins, sbins], weights=rel_power)

        # transform to gradient
        baz_edges = baz_edges / 180 * np.pi

        try:
            fig = self.pylab(get='figure')
        except TypeError:
            fig = plt.figure()
        # add polar and colorbar axes
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        plt.rc('grid',linewidth=0)

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            bars = ax.bar(left=(pi / 2 - (i + 1) * dw) * np.ones(N),
                          height=dh * np.ones(N),
                          width=dw, bottom=dh * np.arange(N),
                          color=cmap(row / hist.max()))
                            

        ax.set_xticks([pi / 2, 0, 3. / 2 * pi, pi])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

        # set slowness limits
        ax.set_ylim(0., self.smax)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=hist.min(), vmax=hist.max()))
        
        plt.show()
        
        print 'Center of Array at latitude %s and longitude %s'%(center_lat, center_lon)

def __snufflings__():
    return [ FK() ]
        

from pyrocko.snuffling import Param, Snuffling, Choice

import numpy as num


def p2o_trace(ptrace, station):
    '''Convert Pyrocko trace to ObsPy trace.'''
    from obspy.core import UTCDateTime
    from obspy.core import Trace as oTrace

    otr = oTrace(
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
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align='center'>FK ANALYSIS</h1>
    <h2 align='center'>based in ObsPy</h2>
    <p>
    <u>Prerequisites</u><br>
    This snuffling requires the ObsPy package which can be found at
        <a href='https://github.com/obspy/obspy/wiki'> ObsPy's github page </a>
    </p>
    <p>
    On way to install this package is to do:
        <PRE>
        git clone git://github.com/obspy/obspy.git
        cd obspy
        sudo python setup.py install
        </PRE>
    <p>
    <u>Usage</u><br>

     - Load station information at startup <br>

     - Zoom into the data until you see only data you desire to analyse or
     use extended markers to selected time regions for analysis<br>

     - Press the 'Run' button <br>

    </p>
    <p>
    The slowness is given in s/km.
    The assumed location of the geometrical center is printed to the terminal.
    </p>
    <p>
    Further information can be gathered from
    <a href="http://docs.obspy.org/master/tutorial/code_snippets/beamforming_fk_analysis.html">
    ObsPy's FK tutorial</a>.
    </p>
    </body>
    </html>
    '''

    def setup(self):
        self.set_name('FK Analysis')
        self.add_parameter(Param('Slowness range[+-]', 'smax', 0.2, 0., 1.))
        self.add_parameter(Param(
            'Number of slowness divisions', 'divisor', 20, 10, 50))
        self.add_parameter(Param(
            'Number of radial sections', 'numberOfFraction', 32, 4, 50))
        self.add_parameter(Param(
            'Length of Sliding Window [s]', 'window_lenth', 1., 0.5, 10.))
        self.add_parameter(Param(
            'Step fraction of Sliding Window [s]','win_frac', 0.05, 0., 10.))
        self.add_parameter(Choice(
            'If sampling rates differ', 'downresample', 'resample',
            ['resample', 'downsample', 'downsample to "target dt"']))
        self.add_parameter(Param('target dt', 'target_dt', 0.2, 0., 10))

        #self.add_parameter(Choice('Units: ','unit','[s/km]',('[s/km]','[s/deg]')))
        self.set_live_update(False)

    def call(self):
        try:
            from obspy.core import UTCDateTime, stream
            from obspy.signal import array_analysis
            from obspy.imaging.cm import obspy_sequential as cmap
        except ImportError as _import_error:
            self.fail('ImportError:\n%s'% _import_error)

        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize
        import matplotlib.dates as mdates
        self.cleanup()
        viewer = self.get_viewer()

        if viewer.lowpass is None or viewer.highpass is None:
            self.fail('highpass and lowpass in viewer must be set!')

        traces = []
        for trs in self.chopper_selected_traces(fallback=True):
            for tr in trs:
                tr.lowpass(2, viewer.lowpass)
                tr.highpass(2, viewer.highpass)

            traces.extend(trs)

        if not traces:
            self.fail('no traces selected')

        if self.downresample == 'resample':
            dt_want = min([t.deltat for t in traces])
            for t in traces:
                t.resample(dt_want)

        elif self.downresample == 'downsample':
            dt_want = max([t.deltat for t in traces])
            for t in traces:
                t.downsample_to(dt_want)

        elif self.downresample == 'downsample to "target dt"':
            for t in traces:
                t.downsample_to(float(self.target_dt))

        tmin = max([t.tmin for t in traces])
        tmax = min([t.tmax for t in traces])
        try:
            obspy_traces = [p2o_trace(
                tr, viewer.get_station(viewer.station_key(tr)))
                            for tr in traces]

        except KeyError:
            self.fail('station information missing')

        st = stream.Stream(traces=obspy_traces)
        center = array_analysis.get_geometry(st, return_center=True)
        center_lon, center_lat, center_ele = center[len(center)-1]

        # Execute sonic
        kwargs = dict(
            sll_x=-self.smax, slm_x=self.smax, sll_y=-self.smax,
            slm_y=self.smax, sl_s=self.smax/self.divisor,
            win_len=self.window_lenth, win_frac=self.win_frac,
            frqlow=viewer.highpass, frqhigh=viewer.lowpass, prewhiten=0,
            semb_thres=-1.0e9, vel_thres=-1.0e9, verbose=True,
            timestamp='mlabday', stime=UTCDateTime(tmin), etime=UTCDateTime(tmax)
        )

        try:
            out = array_analysis.array_processing(st, **kwargs)
        except AttributeError:
            from obspy.signal.array_analysis import sonic
            out = sonic(st, **kwargs)

        pi = num.pi

        # make output human readable, adjust backazimuth to values between 0
        # and 360
        t, rel_power, abs_power, baz, slow = out.T
        baz[baz < 0.0] += 360.

        # choose number of fractions in plot (desirably 360 degree/N is an
        # integer!)
        N = int(self.numberOfFraction)
        abins = num.arange(N + 1) * 360. / N
        sbins = num.linspace(0., self.smax, N + 1)

        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = num.histogram2d(baz, slow,
                bins=[abins, sbins], weights=rel_power)

        # transform to gradient
        baz_edges = baz_edges / 180. * pi

        fig = self.pylab(get='figure')
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        ax.grid(False)

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            bars = ax.bar(left=(pi / 2 - (i + 1) * dw) * num.ones(N),
                          height=dh * num.ones(N),
                          width=dw, bottom=dh * num.arange(N),
                          color=cmap(row / hist.max()))

        ax.set_xticks([pi / 2, 0, 3. / 2 * pi, pi])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])
        ax.set_ylim(0., self.smax)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=hist.min(), vmax=hist.max()))

        fig2 = self.pylab(get='figure')
        labels = ['rel.power', 'abs.power', 'baz', 'slow']
        xlocator = mdates.AutoDateLocator()
        ax = None
        for i, lab in enumerate(labels):
            ax = fig2.add_subplot(4, 1, i + 1, sharex=ax)
            ax.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
                       edgecolors='none', cmap=cmap)
            ax.set_ylabel(lab)
            ax.set_xlim(out[0, 0], out[-1, 0])
            ax.set_ylim(out[:, i + 1].min(), out[:, i + 1].max())
            ax.xaxis.set_tick_params(which='both', direction='in')
            ax.xaxis.set_major_locator(xlocator)
            ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(xlocator))
            if i != 3:
                ax.set_xticklabels([])
        fig2.subplots_adjust(hspace=0.)
        fig2.canvas.draw()
        fig.canvas.draw()

        print 'Center of Array at latitude %s and longitude %s'%(center_lat, center_lon)

def __snufflings__():
    return [ FK() ]

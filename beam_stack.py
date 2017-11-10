from __future__ import print_function

from pyrocko.gui.snuffling import Snuffling, Param, Switch, Choice
from pyrocko.model import Station, dump_stations
from pyrocko import orthodrome as ortho
from pyrocko import util, io, trace
import numpy as num
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from collections import defaultdict
from matplotlib import cm


r_earth = 6371000.785
torad = num.pi/180.
onedeg = r_earth*torad

def to_cartesian(items, reflatlon):
    res = defaultdict()
    for i, item in enumerate(items):

        y, x = ortho.latlon_to_ne(reflatlon, item)
        depth = item.depth
        elevation = item.elevation
        dz = elevation - depth
        lat = item.lat/180.*num.pi
        z = r_earth+dz*num.sin(lat)
        res[item] = (x, y, z)
    return res


class BeamForming(Snuffling):
    '''
    <html>
    <body>
    <h1>Beam Forming</h1>
    Mark a time frame using extended markers. If no time window is selected,
    all traces currently visible in the viewer will be stacked.<br>
    Set a backzimuth and horizontal slowness of the passing wave field.<br>
    If no reference point is defined by 'Center lat' and 'Center lon'
    sliders the geographical center is calculated by taking the average of
    latitudes and longitudes. This can also be done by pressing <b>Set center by
    mean lat/lon</b><br>
    If sampling rates differ and <b>Tread different dt by</b> is set to <b>downsample</b>,
    traces with highest sampling rates (smalles deltat) will be stacked, first. 
    The current stacked trace will then be downsampled to match the next lower sampling
    rate. The stacking will proceed until another lower sampling rate is found, etc.<br>
    Else, if <b>Tread different dt by</b> is set to <b>oversample</b> all traces will by
    resampled in the frequency domain to match the highest occurring sampling rate. 
    <p>
    <b>Tread different dt by </b> - oversample:|downsample.<br>
    <b>pre-normalize by std</b> - normalize traces using their standard deviation.<br>
    <b>multiply 1/[no. of traces]</b> - stacked trace's will be normalized by 
    the number of summed traces.<br>
    <b>Add Shifted Trace</b> - Add time shifted traces to the viewer.<br>
    <b>Plot</b> - map station distribution together with the applied time shifts. The
    arrow indicates the applied back azimuth. Grey dots indicate stations which have not
    been considered in the stacking.</br>
    <b>Save Traces</b> - write stacked traces to mseed file.
    </p>
    </body>
    </html>
    '''
    def setup(self):
        self.set_name("Beam Forming")
        self.add_parameter(Param('Center lat', 'lat_c', 90., -90., 90.,
                                 high_is_none=True))
        self.add_parameter(Param('Center lon', 'lon_c', 180., -180., 180.,
                                 high_is_none=True))
        self.add_parameter(Param('Back azimuth', 'bazi', 0., 0., 360.))
        self.add_parameter(Param('slowness', 'slow', 0.1, 0., 1.))
        self.add_parameter(Choice('slowness unit', 'unit', 's/km',['s/km',
                                                                    's/deg']))
        self.add_parameter(Choice('Treat different dt by', 'diff_dt_treat',
                                  'oversample',['oversample', 'downsample']))
        self.add_parameter(Switch('pre-normalize by std ', 'normalize_std', False))
        self.add_parameter(Switch('multiply 1/[no. of traces]', 'post_normalize', False))
        self.add_parameter(Switch('Add Shifted Traces', 'add_shifted', False))
        self.add_trigger('plot', self.plot)
        self.add_trigger('Save Station', self.save_station)
        self.add_trigger('Save Traces', self.save)
        self.add_trigger('Set center by mean lat/lon', self.set_center_latlon)
        self.station_c = None
        self.z_c = None
        self.stacked_traces = None

    def panel_visibility_changed(self, bool):
        if bool:
            viewer = self.get_viewer()
            self._param_controls['unit'].choosen.connect(
                self.set_slowness_ranges)

    def set_slowness_ranges(self, ident, state):
        if state == 's/km':
            self.set_parameter_range('slow', 0, 1)
            self.set_parameter('slow', self.slow/onedeg*1000.)
        elif state == 's/deg':
            self.set_parameter_range('slow', 0, 100)
            self.set_parameter('slow', self.slow/1000.*onedeg)

    def call(self):
        self.cleanup()
        c_station_id = ('_', 'STK')
        if self.unit == 's/deg':
            slow_factor = 1./onedeg
        elif self.unit == 's/km':
            slow_factor = 1./1000.

        slow = self.slow*slow_factor
        if self.stacked_traces is not None:
            self.add_traces(self.stacked_traces)
        viewer = self.get_viewer()
        if self.station_c:
            viewer.stations.pop(c_station_id)

        stations = self.get_stations()
        if len(stations) == 0:
            self.fail('No station meta information found')

        traces = list(self.chopper_selected_traces(fallback=True))
        traces = [tr for trs in traces for tr in trs ]
        visible_nslcs = [tr.nslc_id for tr in traces]
        stations = [x for x in stations if util.match_nslcs(
            "%s.%s.%s.*" % x.nsl(), visible_nslcs)]
        if not self.lat_c or not self.lon_c or not self.z_c:
            self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(stations)
            self.set_parameter('lat_c', self.lat_c)
            self.set_parameter('lon_c', self.lon_c)

        self.station_c = Station(lat=float(self.lat_c),
                                 lon=float(self.lon_c),
                                 elevation=float(self.z_c),
                                 depth=0.,
                                 name='Array Center',
                                 network=c_station_id[0],
                                 station=c_station_id[1])

        viewer.add_stations([self.station_c])
        lat0 = num.array([self.lat_c]*len(stations))
        lon0 = num.array([self.lon_c]*len(stations))
        lats = num.array([s.lat for s in stations])
        lons = num.array([s.lon for s in stations])
        ns, es = ortho.latlon_to_ne_numpy(lat0, lon0, lats, lons)
        theta = num.float(self.bazi*num.pi/180.)
        R = num.array([[num.cos(theta), -num.sin(theta)],
                        [num.sin(theta), num.cos(theta)]])
        distances = R.dot(num.vstack((es, ns)))[1]
        channels = set()
        self.stacked = {}
        num_stacked = {}
        self.t_shifts = {}
        shifted_traces = []
        taperer = trace.CosFader(xfrac=0.05)
        if self.diff_dt_treat=='downsample':
            traces.sort(key=lambda x: x.deltat)
        elif self.diff_dt_treat=='oversample':
            dts = [t.deltat for t in traces]
            for tr in traces:
                tr.resample(min(dts))

        for tr in traces:
            if tr.nslc_id[:2] == c_station_id:
                continue
            tr = tr.copy(data=True)
            tr.ydata = tr.ydata.astype(num.float64)
            tr.ydata -= tr.ydata.mean(dtype=num.float64)
            tr.taper(taperer)
            try:
                stack_trace = self.stacked[tr.channel]
                num_stacked[tr.channel] += 1
            except KeyError:
                stack_trace = tr.copy(data=True)
                stack_trace.set_ydata(num.zeros(
                    len(stack_trace.get_ydata())))

                stack_trace.set_codes(network=c_station_id[0],
                                      station=c_station_id[1],
                                      location='',
                                      channel=tr.channel)

                self.stacked[tr.channel] = stack_trace
                channels.add(tr.channel)
                num_stacked[tr.channel] = 1

            nslc_id = tr.nslc_id

            try:
                stats = [x for x in stations if util.match_nslc(
                    '%s.%s.%s.*' % x.nsl(), nslc_id)]

                stat = stats[0]
            except IndexError:
                break

            i = stations.index(stat)
            d = distances[i]
            t_shift = d*slow
            tr.shift(t_shift)
            stat = viewer.get_station(tr.nslc_id[:2])
            self.t_shifts[stat] = t_shift
            if self.normalize_std:
                tr.ydata = tr.ydata/tr.ydata.std()

            if num.abs(tr.deltat-stack_trace.deltat)>0.000001:
                if self.diff_dt_treat=='downsample':
                    stack_trace.downsample_to(tr.deltat)
                elif self.diff_dt_treat=='upsample':
                    print('something went wrong with the upsampling, previously')
            stack_trace.add(tr)

            if self.add_shifted:
                tr.set_station('%s_s' % tr.station)
                shifted_traces.append(tr)

        if self.post_normalize:
            for ch, tr in self.stacked.items():
                tr.set_ydata(tr.get_ydata()/num_stacked[ch])

        self.cleanup()

        for ch, tr in self.stacked.items():
            if num_stacked[ch]>1:
                self.add_trace(tr)

        if self.add_shifted:
            self.add_traces(shifted_traces)

    def center_lat_lon(self, stations):
        '''Calculate a mean geographical centre of the array
        using spherical earth'''

        lats = num.zeros(len(stations))
        lons = num.zeros(len(stations))
        elevations = num.zeros(len(stations))
        depths = num.zeros(len(stations))
        for i, s in enumerate(stations):
            lats[i] = s.lat*torad
            lons[i] = s.lon*torad
            depths[i] = s.depth
            elevations[i] = s.elevation

        z = num.mean(elevations-depths)
        return (lats.mean()*180/num.pi, lons.mean()*180/num.pi, z)

    def plot(self):
        stations = self.get_stations()
        res = to_cartesian(stations, self.station_c)
        center_xyz = res[self.station_c]
        x = num.zeros(len(res))
        y = num.zeros(len(res))
        z = num.zeros(len(res))
        sizes = num.zeros(len(res))
        stat_labels = []
        i = 0
        for s, xyz in res.items():
            x[i] = xyz[0]
            y[i] = xyz[1]
            z[i] = xyz[2]

            try:
                sizes[i] = self.t_shifts[s]
                stat_labels.append('%s' % (s.nsl_string()))
            except AttributeError:
                self.fail('Run the snuffling first')
            except KeyError:
                stat_labels.append('%s' % (s.nsl_string()))
                continue
            finally:
                i += 1

        x /= 1000.
        y /= 1000.
        z /= 1000.
        xmax = x.max()
        xmin = x.min()
        ymax = y.max()
        ymin = y.min()

        x_range = num.abs(xmax-xmin)
        y_range = num.abs(ymax-ymin)

        max_range = num.max([x_range, y_range])

        fig = self.pylab(get='figure')
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7])
        ax.set_aspect('equal')
        cmap = cm.get_cmap('bwr')
        ax.scatter(x, y, c=sizes, s=200, cmap=cmap,
                   vmax=num.max(sizes), vmin=-num.max(sizes))
        for i, lab in enumerate(stat_labels):
            ax.text(x[i], y[i], lab, size=14)

        x = x[num.where(sizes==0.)]
        y = y[num.where(sizes==0.)]
        ax.scatter(x, y, c='black', alpha=0.4, s=200)

        ax.arrow(center_xyz[0]/1000.,
                 center_xyz[1]/1000.,
                 -num.sin(self.bazi/180.*num.pi),
                 -num.cos(self.bazi/180.*num.pi),
                 head_width=0.2,
                 head_length=0.2)
        ax.set_ylabel("N-S [km]")
        ax.set_xlabel("E-W [km]")
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=sizes.min(), vmax=sizes.max()))
        fig.canvas.draw()

    def save(self):
        default_fn = 'BeamTraces_baz%s_slow%s.mseed' % (self.bazi, self.slow)
        fn = self.output_filename('Template for output files', default_fn)
        io.save((self.stacked.values()), fn)

    def set_center_latlon(self):
        self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(self.get_stations())
        self.set_parameter('lat_c', self.lat_c)
        self.set_parameter('lon_c', self.lon_c)

    def save_station(self):
        fn = self.output_filename('Save Station')
        dump_stations([self.station_c], fn)

def __snufflings__():
    return [BeamForming()]

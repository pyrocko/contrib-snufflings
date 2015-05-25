from mpl_toolkits.mplot3d import Axes3D
from pyrocko.snuffling import Snuffling, Param, Switch
from pyrocko.model import Station
from pyrocko import orthodrome as ortho
from pyrocko import util, io
import numpy as num
from collections import defaultdict
import matplotlib.pyplot as plt

r_earth = 6371000.785
torad = num.pi/180.
onedeg = r_earth*torad

def to_cartesian(items):
    res = defaultdict()
    latlon00 = ortho.Loc(0.,0.)
    for i, item in enumerate(items):

        y, x = ortho.latlon_to_ne(latlon00, item)
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
    If sampling rates differ, traces with highest sampling rates will be stacked, first. 
    The current stacked trace will then be downsampled to match the next lower sampling
    rate. The stacking will proceed until another lower sampling rate is found, etc.
    <p>
    <b>pre-normalize by std</b> - normalize traces using their standard deviation.<br>
    <b>multiply 1/[no. of traces]</b> - stacked trace's will be normalized by 
    the number of summed traces.<br>
    <b>Add Shifted Trace</b> - Add time shifted traces to the viewer.<br>
    <b>Plot</b> - map station distribution together with the applied time shifts. The
    blue arrow indicates the applied back azimuth.</br>
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
        # min_slow = 0.
        # max_slow = 1.
        # def_slow = 0.1

        self.add_parameter(Param('slowness [s/km]', 'slow', 0.1, 0., 1.))
        # self.add_parameter(Param('slowness [s/deg]',
        #                          'slow_deg',
        #                          def_slow*onedeg,
        #                          min_slow*onedeg,
        #                          max_slow*onedeg))

        self.add_parameter(Switch('pre-normalize by std ', 'normalize_std', False))
        self.add_parameter(Switch('multiply 1/[no. of traces]', 'post_normalize', False))
        self.add_parameter(Switch('Add Shifted Traces', 'add_shifted', False))
        self.add_trigger('plot', self.plot)
        self.add_trigger('Save Traces', self.save)
        self.add_trigger('Set center by mean lat/lon', self.set_center_latlon)
        self.station_c = None
        self.z_c = None
        self.stacked_traces = None


    def call(self):

        self.cleanup()
        if self.stacked_traces is not None:
            self.add_traces(self.stacked_traces)
        viewer = self.get_viewer()
        if self.station_c:
            viewer.stations.pop(('', 'STK'))

        stations = self.get_stations()

        if not self.lat_c or not self.lon_c or not self.z_c:
            self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(stations)
            self.set_parameter('lat_c', self.lat_c)
            self.set_parameter('lon_c', self.lon_c)

        self.station_c = Station(lat=float(self.lat_c),
                                 lon=float(self.lon_c),
                                 elevation=float(self.z_c),
                                 depth=0.,
                                 name='Array Center',
                                 network='',
                                 station='STK')

        viewer.add_stations([self.station_c])

        distances = [ortho.distance_accurate50m(self.station_c, s) for s in
                     stations]

        azirad = self.bazi*torad
        azis = num.array([ortho.azimuth(s, self.station_c) for s in stations])
        azis %= 360.
        azis = azis*torad

        gammas = azis - azirad
        gammas = gammas % (2*num.pi)
        channels = set()
        self.stacked = {}
        num_stacked = {}
        self.t_shifts = {}
        shifted_traces = []
        traces = list(self.chopper_selected_traces(fallback=True))
        traces = [tr for trs in traces for tr in trs ]
        traces.sort(key=lambda x: x.deltat)
        for tr in traces:
            if tr.nslc_id[:3] == ('_', 'STK', ''):
                continue
            tr = tr.copy(data=True)
            tr.ydata = tr.ydata.astype(num.float64)
            tr.ydata -= tr.ydata.mean(dtype=num.float64)
            try:
                stack_trace = self.stacked[tr.channel]
                num_stacked[tr.channel] += 1
            except KeyError:
                stack_trace = tr.copy(data=True)
                stack_trace.set_ydata(num.zeros(
                    len(stack_trace.get_ydata())))

                stack_trace.set_codes(network='_',
                                      station='STK',
                                      location='',
                                      channel=tr.channel)

                self.stacked[tr.channel] = stack_trace
                channels.add(tr.channel)
                num_stacked[tr.channel] = 1

            nslc_id = tr.nslc_id

            try:
                stats = filter(lambda x: util.match_nslc(
                    '%s.%s.%s.*' % x.nsl(), nslc_id), stations)

                stat = stats[0]
            except IndexError:
                break

            i = stations.index(stat)
            gamma = gammas[i]

            d = num.cos(gamma)*distances[i]
            t_shift = d*self.slow/1000.
            tr.shift(-t_shift)
            stat = viewer.get_station(tr.nslc_id[:2])
            self.t_shifts[stat] = -t_shift
            if self.normalize_std:
                tr.ydata = tr.ydata/tr.ydata.std()

            if num.abs(tr.deltat-stack_trace.deltat)>0.000001:
                stack_trace.downsample_to(tr.deltat)
            stack_trace.add(tr)

            if self.add_shifted:
                tr.set_station('%s_s' % tr.station)
                shifted_traces.append(tr)

        #normalize by number of stacked traces:
        if self.post_normalize:
            for ch, tr in self.stacked.items():
                tr.set_ydata(tr.get_ydata()/num_stacked[ch])
        self.stacked_traces = self.stacked.values()
        self.cleanup()
        self.add_traces(self.stacked_traces)
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
        res = to_cartesian(stations)
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
                stat_labels.append('%s\n%1.2f' % (s.nsl_string(), sizes[i]))
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

        fig = plt.figure()
        ax = self.pylab()
        ax.scatter(x, y, c=sizes, s=200, cmap=plt.cm.get_cmap('bwr'),
                   vmax=num.max(sizes), vmin=-num.max(sizes))
        for i, lab in enumerate(stat_labels):
            ax.text(x[i], y[i], lab, size=14)

        ax.arrow(center_xyz[0]/1000.,
                 center_xyz[1]/1000.,
                 num.sin(self.bazi/180.*num.pi),
                 num.cos(self.bazi/180.*num.pi),
                 head_width=1.4,
                 head_length=4)
        #ax.set_xlim([x.mean()-max_range*0.55, x.mean()+max_range*0.55])
        #ax.set_ylim([y.mean()-max_range*0.55, y.mean()+max_range*0.55])
        ax.set_ylabel("N-S [km]")
        ax.set_xlabel("E-W [km]")

    def save(self):
        default_fn = 'BeamTraces_baz%s_slow%s.mseed' % (self.bazi, self.slow)
        fn = self.output_filename('Template for output files', default_fn)
        io.save(self.stacked.values(), fn)

    def set_center_latlon(self):
        self.lat_c, self.lon_c, self.z_c = self.center_lat_lon(self.get_stations())
        self.set_parameter('lat_c', self.lat_c)
        self.set_parameter('lon_c', self.lon_c)

def __snufflings__():
    return [BeamForming()]

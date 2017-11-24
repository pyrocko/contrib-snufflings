from __future__ import print_function

import os
import tempfile
import math
import glob
import numpy as num
import shutil

from pyrocko.gui.snuffling import Snuffling, Switch, Param, Choice
from pyrocko.gui.pile_viewer import EventMarker, PhaseMarker
from pyrocko import util, model, orthodrome
from pyrocko.plot import gmtpy
from subprocess import Popen, PIPE, check_call

deg2rad = math.pi/180.
pjoin = os.path.join


class MPlot(gmtpy.MPlot):

    def pre_draw(self, gmt, widget, scaler):
        gmt.pscoast(D='f', G=(230, 230, 200), W='1p,black', *(widget.JXY() +
                                                              scaler.RB()))


hypo_param_tmpl = '''  hyposat-parameter

*GLOBAL MODEL 2                     : iasp91
GLOBAL MODEL                       : %(global_model)s

LOCAL OR REGIONAL MODEL            : %(locreg_model)s
PHASE INDEX FOR LOCAL MODEL        : 0000
CRUST 5.1                          : %(crust_51)i
CRUST 5.1 PATH                     : ./

OUTPUT OF REGIONAL MODEL (DEF 0)   : 1

STATION FILE                       : stations.dat
STATION CORRECTION FILE            : _

P-VELOCITY TO CORRECT ELEVATION    : %(vp_to_correct_elevation)g
S-VELOCITY TO CORRECT ELEVATION    : %(vs_to_correct_elevation)g

RG GROUP-VELOCITY (DEF 2.5  [km/s]): %(rg_group_velocity)g
LG GROUP-VELOCITY (DEF 3.5  [km/s]): 3.5752

LQ GROUP-VELOCITY (DEF 4.4  [km/s]):  4.4
LR GROUP-VELOCITY (DEF 3.95 [km/s]): 2.85

STARTING SOURCE TIME (EPOCHAL TIME): 0.
STARTING TIME ERROR       [s]      : 0.

STARTING SOURCE DEPTH     [km]     : %(starting_source_depth_km)g
STARTING DEPTH ERROR      [km]     : 50.
DEPTH FLAG (f,b,d,F,B,D)           : b

STARTING SOURCE LATITUDE  [deg]    : 999.
STARTING LATITUDE ERROR   [deg]    : 10.

STARTING SOURCE LONGITUDE [deg]    : 999.
STARTING LONGITUDE ERROR  [deg]    : 10.

MAGNITUDE CALCULATION (DEF 0)      : 1
P-ATTENUATION MODEL (G-R or V-C)   : V-C
S-ATTENUATION MODEL (IASPEI or R-P): R-P

MAXIMUM # OF ITERATIONS            : 80
# TO SEARCH OSCILLATIONS (DEF 4)   : 6

LOCATION ACCURACY [km] (DEFAULT 1) : 1.
CONSTRAIN SOLUTION (0/1)           : 1

CONFIDENCE LEVEL  (68.3 - 99.99 %%) : 95.
EPICENTER ERROR ELLIPSE (DEF 1)    : 1

MAXIMUM AZIMUTH ERROR     [deg]    : 30.
MAXIMUM SLOWNESS ERROR    [s/deg]  : 5.

SLOWNESS [S/DEG] ( 0 = APP. VEL)   : 0

FLAG USING TRAVEL-TIME DIFFERENCES : 1

INPUT FILE NAME (DEF hyposat-in)   : _

OUTPUT FILE NAME (DEF hyposat-out) : _
OUTPUT SWITCH  (YES = 1, DEFAULT)  : 1
OUTPUT LEVEL                       : 4
'''


def nsl_str(nsl):
    return '.'.join(nsl)


def to_min_sec(lat, lon):

    if lat >= 0:
        ns = 'N'
    else:
        ns = 'S'
    if lon >= 0:
        ew = 'E'
    else:
        ew = 'W'

    lat = abs(lat)
    lon = abs(lon)
    dlat = lat - math.floor(lat)
    dlon = lon - math.floor(lon)
    mlat = math.floor(dlat*60)
    mlon = math.floor(dlon*60)
    slat = (dlat*60-mlat)*60
    slon = (dlon*60-mlon)*60

    return '%2i%02i%04.1f%s%3i%02i%04.1f%s' % (math.floor(lat), mlat, slat, ns,
                                               math.floor(lon), mlon, slon, ew)


def ellipse(major, minor, azimuth):
    azimuth = deg2rad*azimuth
    rot = num.array([[math.cos(azimuth), math.sin(azimuth)],
                     [-math.sin(azimuth), math.cos(azimuth)]], dtype=num.float)
    n = 361
    phi = num.linspace(0., 2.*math.pi, n)
    points = num.zeros((n, 2), dtype=num.float)
    points[:, 0] = num.cos(phi)*major
    points[:, 1] = num.sin(phi)*minor
    return num.dot(rot, points.T).T


def ellipse_lat_lon(major, minor, azimuth, lat, lon):
    points = ellipse(major, minor, azimuth)
    return orthodrome.ne_to_latlon(lat, lon, points[:, 0], points[:, 1])


class Hyposat(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Locate seismic sources using HYPOSAT</h1>
    HYPOSAT is a travel time based location software for seismic sources
    developed by Johannes Schweitzer. Download the source codes from
    <a href="ftp://ftp.norsar.no/pub/outgoing/johannes/hyposat/">NORSAR</a> and
    follow the installation instructions provided there.

    <h2>Usage</h2>
    &middot; Invoke snuffler with station meta data.
    <br>
    &middot; Pick <i>P</i> and <i>S</i> phases and define them as
    such by pressing <i>F1</i> and <i>F2</i>, respectively.
    <br>
    &middot; Select markers you want to use in the inversion (e.g.: press
    <i>a</i> to select all markers in the visible time range)
    <br>
    &middot; Press <i>Run</i>
    <p>
    The ouput of HYPOSAT is printed to the terminal. For further information on
    how to tune HYPOSAT read the manual available at
    <a href="ftp://ftp.norsar.no/pub/outgoing/johannes/hyposat/">NORSAR</a>.

    <h2>Reference</h2>
<br>
<i>
    Schweitzer, J. (1997):
    <a href=
        "ftp://ftp.norsar.no/pub/outgoing/johannes/hyposat/schweitzer_1997.pdf">
    HYPOSAT - a new routine to locate seismic events</a>,
   NORSAR Scientific Report 1-97/98, 94-102, NORSAR, Kjeller, Norway,
   November 1997.
   </i>
    </body>
    </html>
    '''
    def setup(self):

        self.hyposat_data_dir = pjoin(self._path, 'hyposat/data')

        locreg_models = [os.path.basename(x) for x in glob.glob(
            pjoin(self.hyposat_data_dir, 'model_*.dat'))]

        crust_51_keys = [
            'Off',
            'For station corrections',
            'For local/regional model',
            'For station corrections, local/regional model and surface\
                        reflection corrections']

        self.crust_51_choices = dict(
            [(b, a) for (a, b) in enumerate(crust_51_keys)])

        self.set_name('HYPOSAT')
        self.add_parameter(
            Choice('Global Model', 'global_model', 'ak135',
                   ('ak135', 'iasp91', 'iasp91a', 'prem', 'jb', 'sp6')))
        self.add_parameter(
            Choice('Local or Regional Model', 'locreg_model', '_', ('_', ) +
                   tuple(locreg_models)))
        self.add_parameter(
            Choice('Use CRUST 5.1', 'crust_51', 'Off', crust_51_keys))
        self.add_trigger('Save', self.save_last_run)
        self.add_parameter(
            Param('P std. deviation [s]', 'p_stddev', 0.1, 0.001, 9.9))
        self.add_parameter(
            Param('S std. deviation [s]', 's_stddev', 0.2, 0.001, 9.9))
        self.add_parameter(
            Param('P-wave velocity to correct elevation',
                  'vp_to_correct_elevation', 3.8, 0., 10.))
        self.add_parameter(
            Param('S-wave velocity to correct elevation',
                  'vs_to_correct_elevation', 2.1, 0., 10.))
        self.add_parameter(
            Param('Starting source depth [km]', 'starting_source_depth_km',
                  10., 0., 600.))
        self.add_parameter(
            Param('Zero level shift [km]', 'zero_level_km',  0., 0., 10.))
        self.add_parameter(
            Param('RG group velocity', 'rg_group_velocity',  2.6, 1., 10.))
        self.add_parameter(
            Switch('Show location plot', 'show_location_plot', False))
        self.set_live_update(False)
        self.pdf_viewer = 'evince'
        self.dir = None

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()

        viewer = self.get_viewer()
        markers = self.get_selected_markers()
        if len(markers) == 0:
            self.fail('No markers selected.')

        event = viewer.get_active_event()

        if len(viewer.stations) == 0:
            self.fail('No station information available.')

        station_phase_to_nslc = {}

        hypo_in = []
        for marker in markers:
            if isinstance(marker, PhaseMarker):
                if marker.get_event() == event:
                    phasename = marker.get_phasename()
                    nslcs = list(marker.nslc_ids)
                    station = nslcs[0][1]

                    station_phase_to_nslc[station, phasename] = nslcs[0]

                    backazi = -1.
                    backazi_stddev = 0.0
                    slowness = -1.
                    slowness_stddev = 0.0
                    period = 0.0
                    amplitude = 0.0
                    flags = 'T__DR_'
                    t = marker.tmin

                    date_str = util.time_to_str(t, '%Y %m %d %H %M %S.3FRAC')
                    if phasename == 'P':
                        t_stddev = self.p_stddev
                    elif phasename == 'S':
                        t_stddev = self.s_stddev

                    hypo_in.append((station, phasename, date_str, t_stddev,
                                    backazi, backazi_stddev, slowness,
                                    slowness_stddev, flags, period, amplitude))

        hypo_in.sort()

        self.dir = tempfile.mkdtemp(prefix='hyposat-%s-' % os.environ['USER'])
        print()
        print('=== Running HYPOSAT ' + '=' * 80)
        print('temp dir: %s' % self.dir)

        fn = pjoin(self.dir, 'hyposat-in')
        f = open(fn, 'w')
        f.write('\n')
        for vals in hypo_in:
            s = '%-5s %-8s %s %5.3f %6.2f %5.2f %5.2f %5.2f %-6s %6.3f %12.2f' % vals
            if len(s) != 96:
                self.fail('Cannot generate input file for HYPOSAT')
                return
            f.write(s+'\n')

        f.close()

        fn = pjoin(self.dir, 'stations.dat')
        f = open(fn, 'w')
        sta_lat_lon = []
        for sta in viewer.stations.values():
            s = '%-5s%1s%s%7.1f' % (
                sta.station, ' ', to_min_sec(sta.lat, sta.lon),
                sta.elevation - self.zero_level_km*1000.)

            f.write(s+'\n')
            sta_lat_lon.append( (sta.lat, sta.lon) )

        f.close()

        params = {}
        params['starting_source_depth_km'] = self.starting_source_depth_km + self.zero_level_km
        params['global_model'] = self.global_model
        params['locreg_model'] = self.locreg_model
        params['vp_to_correct_elevation'] = self.vp_to_correct_elevation
        params['vs_to_correct_elevation'] = self.vs_to_correct_elevation
        params['crust_51'] = self.crust_51_choices[self.crust_51]
        params['rg_group_velocity'] = self.rg_group_velocity

        fn = pjoin(self.dir, 'hyposat-parameter')
        f = open(fn, 'w')
        f.write(hypo_param_tmpl % params)
        f.close()

        old_wd = os.getcwd()
        os.chdir(self.dir)

        env = dict(os.environ)
        env['HYPOSAT_DATA'] = self.hyposat_data_dir

        try:
            p = Popen([pjoin(self._path, 'hyposat/bin/hyposat')], env=env, stdout=PIPE)
        except OSError as e:
            try:
                # try included, compiled version:
                abs_path = os.path.dirname(os.path.abspath(__file__))
                executable = pjoin(abs_path, 'hyposat', 'bin_l', 'hyposat')
                print('abspath', executable)
                p = Popen([executable], env=env, stdout=PIPE)
            except OSError as e:
                import errno
                if e.errno == errno.ENOENT:
                    self.fail('hyposat not found')
                    return
                else:
                    raise e

        (out, err) = p.communicate()
        os.chdir(old_wd)

        fn = pjoin(self.dir, 'hyposat-out')
        f = open(fn, 'r')
        hypo_out = f.read()
        f.close()

        print('HYPOSAT output:\n\n %s' % hypo_out)
        evhead = 'T0 LAT LON Z VPVS DLAT DLON DZ DT0 DVPVS DEF RMS'.split()
        phhead = 'Stat Delta Azi Phase [used] Onset time Res Baz Res Rayp Res Used'.split()
        ellipsehead = 'Epicenter error ellipse:'.split()
        state = 0
        source_time = None
        phmarks = []
        kind = 1
        event_markers = []
        ellipse_major = None
        ellipse_minor = None
        ellipse_azimuth = None
        for line in hypo_out.splitlines():
            if state == 0:
                if line.split() == evhead:
                    state = 1
                elif line.split() == phhead:
                    state = 2
                elif line.split() == ellipsehead:
                    state = 3
                elif line.lstrip().startswith('Source time  :'):
                    toks = line.split()
                    datestr = ' '.join(toks[3:9])
                    source_date_str = ' '.join(toks[3:6])
                    source_time = util.str_to_time(datestr, format='%Y %m %d %H %M %S.3FRAC')

            elif state == 1:
                toks = line.split()
                datestr = ' '.join(toks[:4])
                t = util.str_to_time(datestr, format='%Y-%m-%d %H %M %S.3FRAC')
                lat = float(toks[4])
                lon = float(toks[5])
                depth = float(toks[6])*1000.
                event = model.Event(lat, lon, t, depth=depth-(self.zero_level_km*1000.), name='HYPOSAT-%i' % kind)
                event.ellipse_major = ellipse_major
                event.ellipse_minor = ellipse_minor
                event.ellipse_azimuth = ellipse_azimuth
                evmark = EventMarker(event)
                event_markers.append(evmark)
                evmark.set_kind(kind)
                self.add_marker(evmark)

                for phmark in phmarks:
                    phmark.set_kind(kind)
                    phmark.set_event(event)
                    self.add_marker(phmark)

                phmarks = []
                kind += 1
                state = 0

            elif state == 2:
                toks = line[:58].split()
                if len(toks) >= 8:
                    have_used = 0
                    if len(toks) == 9:
                        have_used = 1
                    station = toks[0]
                    phase = toks[3]
                    if have_used:
                        used_phase = toks[4]
                    else:
                        used_phase = phase
                    residual = float(toks[7+have_used])
                    timestr = ' '.join(toks[4+have_used:7+have_used])
                    datestr = source_date_str + ' ' + timestr
                    t = util.str_to_time(datestr, format='%Y %m %d %H %M %S.3FRAC') - residual
                    if t - source_time < - 12*3600.:
                        t += 24*3600.

                    if t - source_time > 12*3600.:
                        t -= 24*3600.

                    nslc = station_phase_to_nslc[station, phase]
                    phmarks.append(PhaseMarker([nslc], t, t, 0, phasename=used_phase))

                else:
                    if len(toks) == 0 and phmarks:
                        state = 0

            elif state == 3:
                toks = line.split()
                if len(toks) == 0:
                    state = 0
                elif toks[0] == 'Major':
                    ellipse_major = float(toks[2])*1000.
                    ellipse_minor = float(toks[6])*1000.
                elif toks[0] == 'Azimuth:':
                    ellipse_azimuth = float(toks[1])

        print('='*100)
        print()
        if self.show_location_plot:
            cm = gmtpy.cm
            p = MPlot(width=15*cm, height=15*cm, gmtconfig={ 'PLOT_DEGREE_FORMAT':'D' })
            p.set(xspace=0.05, yspace=0.05)
            lats, lons = num.array(sta_lat_lon, dtype=num.float).T

           # p.plot((lons,lats), '-Ggray -W0.5p,black -St10p')
            stations_used = [ x[0] for x in hypo_in ]
            for sta in viewer.stations.values():
                if sta.station in stations_used:
                    p.text([sta.lon, sta.lat, sta.station], justify='CT', offset=(0,-5))
                    p.plot(([sta.lon], [sta.lat]), '-Gblack -W0.5p,black -St10p')

            for evmark in event_markers:
                color = evmark.select_color(evmark.color_b)
                ev = evmark.get_event()
                p.plot([[ev.lon], [ev.lat]], '-Sa20p -G%s' % gmtpy.color(color))
                #p.plot([[ev.lon], [ev.lat]], '-Sa20p -Gred -Wblack')
                elat, elon = ellipse_lat_lon(ev.ellipse_major, ev.ellipse_minor, ev.ellipse_azimuth, ev.lat, ev.lon)
                p.plot([elon, elat], '-W0.5p,%s' % gmtpy.color(color))

            fn_plot = pjoin(self.dir, 'location.png')
            try:
                p.save(fn_plot)
                f = self.pixmap_frame()
                f.load_pixmap(fn_plot)
            except gmtpy.GMTInstallationProblem as e:
                msg = "%s\nFailed to start GMT" % e
                self.fail(msg)

    def save_last_run(self):
        if not self.dir:
            self.fail('Run first')
        outdir = self.output_filename(caption='Save directory')
        shutil.move(self.dir, outdir)


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [Hyposat()]

import numpy as num
import time
from matplotlib.animation import FuncAnimation

from pyrocko.snuffling import Param, Snuffling, Choice, Switch
from scipy.signal import fftconvolve, lfilter, hilbert
from scipy.interpolate import UnivariateSpline
from pyrocko import orthodrome as ortho
from pyrocko import parstack
from pyrocko import util
from pyrocko import trace
from pyrocko import model
import logging


logger = logging.getLogger('pyrocko.snufflings.fk_parstack.py')
d2r = num.pi/180.
km = 1000.


def search_max_block(n_maxsearch, data):
    '''
    Find indices of maxima of *data* in groups of *n_maxsearch*
    samples.

    returns an array of indices.

    If length of *data* is not a multiple of *n_maxsearch*, *data* will be
    padded.
    '''
    n = len(data)
    n_dim2 = (int(n/n_maxsearch)+1) * n_maxsearch
    n_missing = n_dim2 - n
    if n % n_maxsearch != 0:
        a = num.pad(data, [0, n_missing], mode='minimum')
    else:
        a = data

    return num.argmax(a.reshape((-1, n_maxsearch)), axis=1) +\
        num.arange(0, (n+n_missing)/n_maxsearch) * n_maxsearch


def instantaneous_phase(signal):
    analytic_signal = hilbert(signal)
    return num.unwrap(num.angle(analytic_signal))


def get_instantaneous_frequency(signal, fs):
    inst_phase = instantaneous_phase(signal)
    return (num.diff(inst_phase) / (2.0*num.pi) * fs)


def get_center_station(stations, select_closest=False):
    ''' gravitations center of *stations* list.'''
    n = len(stations)
    slats = num.empty(n)
    slons = num.empty(n)
    for i, s in enumerate(stations):
        slats[i] = s.lat
        slons[i] = s.lon

    center_lat, center_lon = ortho.geographic_midpoint(slats, slons)

    if select_closest:
        center_lats = num.ones(n)*center_lat
        center_lons = num.ones(n)*center_lon
        dists = ortho.distance_accurate50m_numpy(
            center_lats, center_lons, slats, slons)

        return stations[num.argmin(dists)]
    else:
        return model.Station(center_lat, center_lon)


def get_theoretical_backazimuth(event, stations, center_station):
    return (ortho.azimuth_numpy(
        event.lat, event.lon, center_station.lat, center_station.lon)
            + 180.) % 360.


def get_shifts(stations, center_station, bazis, slownesses):
    ''' shape = (len(bazi), len(slow))'''
    lats = num.array([s.lat for s in stations])
    lons = num.array([s.lon for s in stations])
    lat0 = num.array([center_station.lat] * len(stations))
    lon0 = num.array([center_station.lon] * len(stations))
    ns, es = ortho.latlon_to_ne_numpy(lat0, lon0, lats, lons)
    station_vector = num.array((ns, es)).T
    bazis = bazis * d2r
    shifts = num.zeros((len(bazis)*len(slownesses), len(stations)))
    ishift = 0
    for ibazi, bazi in enumerate(bazis):
        s_vector = num.array((num.cos(bazi), num.sin(bazi)))
        for islowness, slowness in enumerate(slownesses):
            shifts[ishift] = station_vector.dot(s_vector) * slowness
            ishift += 1

    return shifts


def to_db(d):
    return 10*num.log10(d/num.max(d))


def lowpass_array(ydata_array, deltat, order, corner, demean=True, axis=1):
    '''
    Apply butterworth highpass to the trace.

    :param order: order of the filter
    :param corner: corner frequency of the filter

    Mean is removed before filtering.
    '''

    (b, a) = trace._get_cached_filter_coefs(
        order, [corner*2.0*deltat], btype='low')
    data = ydata_array.astype(num.float64)
    if len(a) != order+1 or len(b) != order+1:
        logger.warn(
            'Erroneous filter coefficients returned by '
            'scipy.signal.butter(). You may need to downsample the '
            'signal before filtering.')
    if demean:
        data -= num.mean(data, axis=1)[None].T

    return lfilter(b, a, data)


def highpass_array(ydata_array, deltat, order, corner, demean=True, axis=1):
    '''
    Apply butterworth highpass to the trace.

    :param order: order of the filter
    :param corner: corner frequency of the filter

    Mean is removed before filtering.
    '''
    (b, a) = trace._get_cached_filter_coefs(
        order, [corner*2.0*deltat], btype='high')
    data = ydata_array.astype(num.float64)
    if len(a) != order+1 or len(b) != order+1:
        logger.warn(
            'Erroneous filter coefficients returned by '
            'scipy.signal.butter(). You may need to downsample the '
            'signal before filtering.')
    if demean:
        data -= num.mean(data, axis=1)[None].T

    return lfilter(b, a, data)


def value_to_index(value, range_min, range_max, range_delta, clip=True):
    ''' map a(n array of) *values* to its' index in a continuous data range
    defined by *range_min*, *range_max* and *range_delta*.
    '''
    indices = num.round((value-range_min)/range_delta)
    if clip:
        indices = num.clip(indices, 0, (range_max-range_min)/range_delta)

    return num.asarray(indices, dtype=num.int)


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
    <p>
    Performs delay and sum in the time domain.
    <u>Usage</u><br>

     - Load station information at startup <br>

     - Zoom into the data until you see only data you desire to analyse or
     use extended markers to selected time regions for analysis<br>

     - Press the 'Run' button <br>

    </p>
    The slowness is in units s/km.
    <p>
    if <b>Show</n> is activated, three images will be genereated for each
    processing block: A polar plot which shows the maximum coherence found
    along the time axis for all slownesses and back-azimuths. The blue dot
    within that figure indicates the position of the maximum. Based on that
    slowness/back-azimuth the other two coherence maps are generated which
    show the coherence in the slowness and back-azimuth domain for that
    specific maximum of that processing block.<br>

    Picinbono, et. al, 1997, On Instantaneous Amplitude and Phase of Signals,
    552 IEEE TRANSACTIONS ON SIGNAL PROCESSING, 45, 3, March 1997
    </body>
    </html>
    '''

    def setup(self):
        self.set_name('FK (parstack)')
        self.add_parameter(Param(
            'max slowness [s/km]', 'slowness_max', 0.2, 0., 1.))
        self.add_parameter(Param(
            'min slowness [s/km]', 'slowness_min', 0.01, 0., 1.))
        self.add_parameter(Param(
            'delta slowness [s/km]', 'slowness_delta', 0.002, 0., 0.2))
        self.add_parameter(Param(
            'delta backazimut', 'delta_bazi', 2, 1, 20))
        self.add_parameter(Param(
            'Increment [s]', 'tinc', 10., 0.5, 10.,
            high_is_none=True))
        self.add_parameter(Param(
            'Smoothing length [N]', 'ntaper', 0, 0, 30, low_is_none=True))
        self.add_parameter(Choice(
            'Use channels', 'want_channel', '*',
            ['*', '*Z', '*E', '*N', 'SHZ', 'BHZ', 'p0']))
        self.add_parameter(
            Choice('method', 'method', 'stack', ['stack', 'correlate'])
        )
        self.add_parameter(Switch('Show', 'want_all', True))
        self.add_parameter(Switch('Phase weighted stack', 'want_pws', False))
        self.add_trigger('Clear Figures', self.cleanup_figures)
        self.set_live_update(False)
        self.irun = 0
        self.figure_frames = []
        self.figs2draw = []

    def cleanup_figures(self):
        '''close all figures.'''
        parent = self._panel_parent
        for fframe in self.figure_frames:
            parent.remove_tab(fframe)
        self.figure_frames = []

    def new_figure(self, title=''):
        '''Return a new Figure instance'''
        fig_frame = self.pylab(name='FK: %s (%i)' %
                          (title, self.irun), get='figure_frame')
        self.figure_frames.append(fig_frame)
        self.figs2draw.append(fig_frame.gcf())
        return self.figs2draw[-1]

    def draw_figures(self):
        ''' Draw all new figures and clear list.'''
        for fig in self.figs2draw:
            fig.canvas.draw()

        self.figs2draw = []

    def call(self):

        self.cleanup()
        figs = []
        azi_theo = None
        method = {'stack': 0,
                  'correlate': 2}[self.method]

        bazis = num.arange(0., 360.+self.delta_bazi, self.delta_bazi)
        slownesses = num.arange(self.slowness_min/km,
                                self.slowness_max/km,
                                self.slowness_delta/km)
        n_bazis = len(bazis)
        n_slow = len(slownesses)

        viewer = self.get_viewer()
        event = viewer.get_active_event()

        stations = self.get_stations()
        stations_dict = dict(zip([viewer.station_key(s) for s in stations],
                                 stations))

        traces_pile = self.get_pile()
        deltats = traces_pile.deltats.keys()
        if len(deltats) > 1:
            self.fail('sampling rates differ in dataset')
        else:
            deltat_cf = deltats[0]

        tinc_use = self.get_tinc_use(precision=deltat_cf)

        if self.ntaper:
            taper = num.hanning(int(self.ntaper))
        else:
            taper = None

        frames = None
        t1 = time.time()

        # make sure that only visible stations are used
        use_stations = stations
        center_station = get_center_station(use_stations, select_closest=True)
        print 'Center station: ', center_station
        shift_table = get_shifts(
            stations=use_stations,
            center_station=center_station,
            bazis=bazis,
            slownesses=slownesses)

        shifts = num.round(shift_table / deltat_cf).astype(num.int32)

        # padding from maximum shift of traces:
        npad = num.max(num.abs(shifts))
        tpad = npad * deltat_cf

        # additional padding for cross over fading
        npad_fade = 0
        tpad_fade = npad_fade * deltat_cf

        npad += npad_fade
        tpad += tpad_fade

        frames = None
        tinc_add = tinc_use or 0

        def trace_selector(x):
            return util.match_nslc('*.*.*.%s' % self.want_channel, x.nslc_id)

        for traces in self.chopper_selected_traces(
                tinc=tinc_use, tpad=tpad, fallback=True,
                want_incomplete=False, trace_selector=trace_selector):

            if len(traces) == 0:
                self.fail('No traces matched')
                continue

            # should be correct
            t_min = traces[0].tmin
            t_max = traces[0].tmax

            use_stations = []
            for tr in traces:
                try:
                    use_stations.append(stations_dict[viewer.station_key(tr)])
                except KeyError:
                    self.fail('no trace %s' % ('.'.join(tr.nslc_id)))

            shift_table = get_shifts(
                stations=use_stations,
                center_station=center_station,
                bazis=bazis,
                slownesses=slownesses)

            shifts = num.round(shift_table / deltat_cf).astype(num.int32)

            wmin = traces[0].tmin
            wmax = wmin + tinc_add

            iwmin = int(round((wmin-wmin) / deltat_cf))
            iwmax = int(round((wmax-wmin) / deltat_cf))
            lengthout = iwmax - iwmin
            arrays = num.zeros((len(traces), lengthout + npad*2))

            for itr, tr in enumerate(traces):
                tr = tr.copy()
                if viewer.highpass:
                    tr.highpass(4, viewer.highpass, demean=True)
                else:
                    tr.ydata = num.asarray(
                        tr.ydata, dtype=num.float) - num.mean(tr.ydata)
                if viewer.lowpass:
                    tr.lowpass(4, viewer.lowpass)

                arrays[itr] = tr.get_ydata()

            # if viewer.highpass:
            #     arrays = highpass_array(
            #            arrays, deltat_cf, 4, viewer.highpass)
            # if viewer.lowpass:
            #     arrays = lowpass_array(arrays, deltat_cf, 4, viewer.lowpass)

            _arrays = []
            for itr, tr in enumerate(traces):
                if taper is not None:
                    ydata = fftconvolve(arrays[itr], taper, mode='same')
                else:
                    ydata = arrays[itr]
                _arrays.append(num.asarray(ydata, dtype=num.float64))
            arrays = _arrays

            offsets = num.array(
                [int(round((tr.tmin-wmin) / deltat_cf)) for tr in traces],
                dtype=num.int32)

            ngridpoints = len(bazis)*len(slownesses)
            weights = num.ones((ngridpoints, len(traces)))

            frames, ioff = parstack.parstack(
                arrays, offsets, shifts, weights, method,
                offsetout=iwmin,
                lengthout=lengthout,
                result=frames,
                impl='openmp')

            # theoretical bazi
            if event is not None:
                azi_theo = get_theoretical_backazimuth(
                    event, use_stations, center_station)
                print('theoretical azimuth %s degrees' % (azi_theo))

            print('processing time: %s seconds' % (time.time()-t1))

            if frames is None:
                self.fail('Could not process data!')
                return

            frames_reshaped = frames.reshape((n_bazis, n_slow, lengthout))
            times = num.linspace(t_min-tpad_fade, t_max+tpad_fade, lengthout)
            max_powers = num.max(frames, axis=0)

            # power maxima in blocks
            i_max_blocked = search_max_block(
                n_maxsearch=npad, data=max_powers)

            if self.want_all:

                # ---------------------------------------------------------
                # maxima search
                # ---------------------------------------------------------
                fig1 = self.new_figure('Max Power')
                nsubplots = 1
                ax = fig1.add_subplot(nsubplots, 1, 1)
                ax.plot(num.max(frames, axis=0))
                _argmax = num.argmax(frames, axis=0)
                imax_bazi_all, imax_slow_all = num.unravel_index(
                    _argmax, dims=(n_bazis, n_slow))

                max_powers += (num.min(max_powers)*-1)
                max_powers /= num.max(max_powers)
                max_powers *= 10.  # Maximum has pixel size of 10
                max_powers *= max_powers
                block_max_times = times[i_max_blocked]
                # --------------------------------------------------------------
                # coherence maps
                # --------------------------------------------------------------

                max_time = num.amax(frames, axis=0)
                imax_time = num.argmax(max_time)
                best_frame = num.amax(frames, axis=1)
                imax_bazi_slow = num.argmax(best_frame)
                imax_bazi, imax_slow = num.unravel_index(
                    num.argmax(best_frame),
                    dims=(n_bazis, n_slow))

                fig2 = self.new_figure('Slowness')
                data = frames_reshaped[imax_bazi, :, :]
                data_max = num.amax(frames_reshaped, axis=0)

                ax = fig2.add_subplot(211)
                ax.set_title('Global maximum slize')
                ax.set_ylabel('slowness [s/km]')
                ax.pcolormesh(times, slownesses*km, data)

                ax = fig2.add_subplot(212, sharex=ax, sharey=ax)
                ax.set_ylabel('slowness [s/km]')
                ax.pcolormesh(times, slownesses*km, data_max)
                ax.set_title('Maximum')

                # highlight block maxima
                local_max_slow = slownesses[imax_slow_all]*km
                ax.plot(block_max_times, local_max_slow[i_max_blocked], 'wo')

                # spline
                spline_slow = UnivariateSpline(
                    block_max_times,
                    local_max_slow[i_max_blocked],
                    w=max_powers[i_max_blocked],
                )

                slow_fitted = spline_slow(times)
                ax.plot(times, num.clip(
                    slow_fitted, self.slowness_min, self.slowness_max)
                )

                fig3 = self.new_figure('Back-Azimuth')
                data = frames_reshaped[:, imax_slow, :]
                data_max = num.amax(frames_reshaped, axis=1)

                ax = fig3.add_subplot(211, sharex=ax)
                ax.set_title('Global maximum slize')
                ax.set_ylabel('back-azimuth')
                ax.pcolormesh(times, bazis, data)
                ax.plot(times[imax_time], bazis[imax_bazi], 'b.')

                ax = fig3.add_subplot(212, sharex=ax, sharey=ax)
                ax.set_ylabel('back-azimuth')
                ax.set_title('Maximum')
                ax.pcolormesh(times, bazis, data_max)

                # highlight block maxima
                local_max_bazi = bazis[imax_bazi_all]
                ax.plot(block_max_times, local_max_bazi[i_max_blocked], 'wo')

                spline_bazi = UnivariateSpline(
                    block_max_times,
                    local_max_bazi[i_max_blocked],
                    w=max_powers[i_max_blocked],
                    k=3,
                    s=4e7
                )

                bazi_fitted = spline_bazi(times)
                ax.plot(times, num.clip(bazi_fitted, 0, 360.))

                i_bazi_fitted = value_to_index(
                    bazi_fitted, 0., 360., self.delta_bazi)
                i_slow_fitted = value_to_index(
                    slow_fitted, self.slowness_min, self.slowness_max,
                    self.slowness_delta)

                i_shift = num.ravel_multi_index(
                    num.vstack((i_bazi_fitted, i_slow_fitted)),
                    (n_bazis, n_slow),
                )

                print 'XX', lengthout
                stack_trace = num.zeros(lengthout)
                i_base = num.arange(lengthout, dtype=num.int)
                for itr, tr in enumerate(traces):
                    print 'XXX', len(tr.ydata)
                    isorting = num.clip(
                        i_base-shifts[i_shift, itr], 0, lengthout)
                    stack_trace += tr.ydata[isorting]

                # xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
                # ax.xaxis.set_major_formatter(xfmt)
                # fig.autofmt_xdate()
                # fig.subplots_adjust(hspace=0)

                semblance = best_frame.reshape((n_bazis, n_slow))

                fig4 = self.new_figure('Max')
                theta, r = num.meshgrid(bazis, slownesses)
                theta *= (num.pi/180.)

                ax = fig4.add_subplot(111, projection='polar')
                m = ax.pcolormesh(theta.T, r.T*km, to_db(semblance))

                ax.plot(bazis[imax_bazi]*d2r, slownesses[imax_slow]*km, 'o')

                bazi_max = bazis[imax_bazi]*d2r
                slow_max = slownesses[imax_slow]*km
                ax.plot(bazi_max, slow_max, 'b.')
                ax.text(0.5, 0.01, 'Maximum at %s degrees, %s s/km' %
                        (num.round(bazi_max, 1), slow_max),
                        transform=fig4.transFigure,
                        horizontalalignment='center',
                        verticalalignment='bottom')

                if azi_theo:
                    ax.arrow(azi_theo/180.*num.pi, num.min(slownesses), 0,
                             num.max(slownesses), alpha=0.5, width=0.015,
                             edgecolor='black', facecolor='green', lw=2,
                             zorder=5)

                self.adjust_polar_axis(ax)
                fig4.colorbar(m)

                # ---------------------------------------------------------
                # CF and beam forming
                # ---------------------------------------------------------
                fig5 = self.new_figure('Beam')
                nsubplots = 4
                nsubplots += self.want_pws

                ax_raw = fig5.add_subplot(nsubplots, 1, 1)
                ax_shifted = fig5.add_subplot(nsubplots, 1, 2)
                ax_beam = fig5.add_subplot(nsubplots, 1, 3)
                ax_beam_new = fig5.add_subplot(nsubplots, 1, 4)

                axkwargs = dict(alpha=0.3, linewidth=0.3, color='grey')

                ybeam = num.zeros(lengthout)
                ybeam_weighted = num.zeros(lengthout)
                for i, (shift, array) in enumerate(zip(shifts.T, arrays)):
                    ax_raw.plot(times, array[npad: -npad], **axkwargs)
                    ishift = shift[imax_bazi_slow]
                    ax_shifted.plot(
                        times, array[npad-ishift: -npad-ishift], **axkwargs)

                    ydata = traces[i].get_ydata()[npad-ishift: -npad-ishift]
                    ybeam += ydata

                    # calculate phase weighting
                    if self.want_pws:
                        ph_inst = instantaneous_phase(ydata)
                        ybeam_weighted += num.abs(num.exp(ph_inst))**4

                ax_beam_new.plot(stack_trace)
                ax_beam_new.set_title('continuous mode')
                # ax_beam.plot(times, ybeam, color='black')
                ax_beam.plot(ybeam, color='black')
                ax_raw.set_title('Characteristic Function')
                ax_shifted.set_title('Shifted CF')
                ax_beam.set_title('Linear Stack')

                if self.want_pws:
                    ax_playground = fig5.add_subplot(nsubplots, 1, 4)
                    ax_playground.plot(ybeam*ybeam_weighted/len(arrays))
                    ax_playground.set_title('Phase Weighted Stack')

                # beam_tr = trace.Trace(
                #     tmin=t_min+tpad, ydata=ybeam, deltat=deltat_cf)
                beam_tr = trace.Trace(
                    tmin=t_min+tpad, ydata=stack_trace, deltat=deltat_cf)

                # -----------------------------------------------------------
                # polar movie:
                # -----------------------------------------------------------
                fig6 = self.new_figure('Beam')
                self.polar_movie(
                    fig=fig6,
                    frames=frames,
                    times=times,
                    theta=theta.T,
                    r=r.T*km,
                    nth_frame=2,
                    n_bazis=n_bazis,
                    n_slow=n_slow,
                )

                self.add_trace(beam_tr)
                self.draw_figures()

                self.irun += 1

    def polar_movie(self, fig, frames, times, theta, r, nth_frame,
                    n_bazis, n_slow):
        frame_artists = []
        # progress_artists = []
        ax = fig.add_subplot(111, projection='polar')
        iframe_min = 0
        iframe_max = len(times)-1
        vmin = num.min(frames)
        vmax = num.max(frames)
        self.adjust_polar_axis(ax)

        def update(iframe):
            # if iframe is not None:
            if False:
                frame = frames[:, iframe]
                '''
                if not progress_artists:
                    progress_artists[:] = [axes2.axvline(
                        tmin_frames - t0 + deltat_cf * iframe,
                        color=scolor('scarletred3'),
                        alpha=0.5,
                        lw=2.)]

                else:
                    progress_artists[0].set_xdata(
                        tmin_frames - t0 + deltat_cf * iframe)
                '''
            else:
                frame = frames[:, iframe].reshape((n_bazis, n_slow))

            frame_artists[:] = [ax.pcolormesh(theta, r, frame, vmin=vmin,
                                              vmax=vmax)]

            # return frame_artists + progress_artists + static_artists
            return frame_artists

        axf = FuncAnimation(  # noqa
            fig, update,
            frames=list(
                xrange(iframe_min, iframe_max+1))[::nth_frame] + [None],
            interval=20.,
            repeat=False,
            blit=True)

        fig.canvas.draw()

    def adjust_polar_axis(self, ax):
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_xticks([0, num.pi/2., num.pi, 3*num.pi/2])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

    def get_tinc_use(self, precision=1.):
        '''
        Set increment time for data processing.
        '''
        if self.tinc is not None:
            tinc = self.tinc
        else:
            tmin, tmax = self.get_selected_time_range(fallback=True)
            if tmin == tmax:
                # selected event marker
                tmin, tmax = self.get_viewer().get_time_range()
            tinc = tmax-tmin
        return num.floor(tinc/precision) * precision


def __snufflings__():
    return [FK()]

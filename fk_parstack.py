import numpy as num
import time

from pyrocko.snuffling import Param, Snuffling, Choice, Switch
from scipy.signal import fftconvolve, argrelextrema, lfilter, hilbert
from pyrocko import orthodrome as ortho
from pyrocko import parstack
from pyrocko import util
from pyrocko import trace
import logging


logger = logging.getLogger('pyrocko.snufflings.fk_parstack.py')
d2r = num.pi/180.
km = 1000.


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
        dists = ortho.distance_accurate50m_numpy(center_lats, center_lons, slats, slons)
        return stations[num.argmin(dists)]
    else:
        return model.Station(center_lat, center_lon)


def get_theoretical_backazimuth(event, stations, center_station):
    return (ortho.azimuth_numpy(event.lat, event.lon, center_station.lat, center_station.lon)\
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

    Picinbono, et. al, 1997, On Instantaneous Amplitude and Phase of Signals, 552
    IEEE TRANSACTIONS ON SIGNAL PROCESSING, 45, 3, March 1997
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
            'delta slowness [s/km]', 'slowness_delta', 0.001, 0., 1.))
        self.add_parameter(Param(
            'delta backazimut', 'delta_bazi', 1, 1, 20))
        self.add_parameter(Param(
            'Increment [s]', 'tinc', 10., 0.5, 10.,
            high_is_none=True))
        self.add_parameter(Param(
            'Smoothing length [N]', 'ntaper', 0, 0, 30, low_is_none=True))
        self.add_parameter(Param('Window length [s]', 'win_length', 2, 0, 30))
        self.add_parameter(Choice(
            'Use channels', 'want_channel', '*',
            ['*', '*Z', '*E', '*N', 'SHZ', 'BHZ', 'p0']))
        self.add_parameter(Switch('Show', 'want_all', True))
        self.add_parameter(Switch('Phase weighted stack', 'want_pws', False))
        self.set_live_update(False)
        self.irun = 0

    def call(self):

        self.cleanup()
        fig1 = None
        azi_theo = None
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

        if self.want_all:
            method = 0
        else:
            method = 1

        def trace_selector(x):
            return util.match_nslc('*.*.*.%s' % self.want_channel, x.nslc_id)

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

            # add lowpass/highpass array!!
            for itr, tr in enumerate(traces):
                tr = tr.copy()
                if viewer.highpass:
                    tr.highpass(4, viewer.highpass, demean=True)
                else:
                    self.fail('Main controls highpass is required')
                if viewer.lowpass:
                    tr.lowpass(4, viewer.lowpass)

                arrays[itr] = tr.get_ydata()

            ntraces = len(traces)

            #if viewer.highpass:
            #    arrays = highpass_array(arrays, deltat_cf, 4, viewer.highpass)
            #if viewer.lowpass:
            #    arrays = lowpass_array(arrays, deltat_cf, 4, viewer.lowpass)
            #arrays *= arrays

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
                azi_theo = get_theoretical_backazimuth(event, use_stations, center_station)
                print('theoretical azimuth %s degrees' % (azi_theo))

            print('processing time: %s seconds' % (time.time()-t1))

            if frames is None:
                self.fail('Could not process data!')
                return

            frames_reshaped = frames.reshape((n_bazis, n_slow, lengthout))
            times = num.linspace(t_min-tpad_fade, t_max+tpad_fade, lengthout)

            # --------------------------------------------------------------
            # Windowed maxima search
            # --------------------------------------------------------------
            #fig1 = self.pylab(name='FK: Power (%i)'%self.irun, get='figure')
            #power_vs_t = num.amax(num.amax(frames, axis=1), axis=0)
            ##print power_vs_t
            #axs = []
            #for i in range(2):
            #    axs.append(fig1.add_subplot(2, 1, i+1))

            #plot_back_slow_time(
                #bazis=bazis, slownesses=slownesses, times=times, frames=frames_reshaped,
            #    bazis=bazis, slownesses=slownesses, times=times, frames=frames,
            #    axs=axs, window_length=self.win_length*deltat_cf)
            #ax.plot(power_vs_t)
            #i_best_t = num.argmax(num.argmax(frames, axis=0), axis=0)
            #i_best_bazi_slow = num.argmax(power_vs_t)
            #imax_bazi, imax_slow = num.unravel_index(
            #i_best_bazi_slow, dims=(n_bazis, n_slow))

            # --------------------------------------------------------------
            # coherence maps
            # --------------------------------------------------------------
            if self.want_all:

                max_time = num.amax(frames, axis=0)
                imax_time = num.argmax(max_time)
                best_frame = num.amax(frames, axis=1)
                imax_bazi_slow = num.argmax(best_frame)
                imax_bazi, imax_slow = num.unravel_index(
                    num.argmax(best_frame),
                    dims=(n_bazis, n_slow))

                fig = self.pylab(name='FK: Slowness (%i)' % self.irun, get='figure')

                data = frames_reshaped[imax_bazi, :, :]
                grad = num.abs(num.gradient(data, edge_order=2)[1])**2

                ax = fig.add_subplot(211)
                ax.set_ylabel('slowness [s/km]')
                ax.pcolormesh(times, slownesses*km, data)
                ax.plot(times[imax_time], slownesses[imax_slow]*km, 'b.')

                ax = fig.add_subplot(212, sharey=ax)
                ax.set_ylabel('abs(grad(slowness))')
                ax.set_xlabel('time')
                ax.pcolormesh(times, slownesses*km, grad)

                fig = self.pylab(name='FK: Back-Azimuth (%i)' % self.irun, get='figure')
                data = frames_reshaped[:, imax_slow, :]
                grad = num.gradient(data, edge_order=2)[1]**2

                ax = fig.add_subplot(211, sharex=ax)
                ax.set_ylabel('back-azimuth')
                ax.pcolormesh(times, bazis, data)
                ax.plot(times[imax_time], bazis[imax_bazi], 'b.')
                ax.set_title('')

                ax = fig.add_subplot(212, sharex=ax, sharey=ax)
                ax.set_ylabel('abs(grad(back-azimuth))')
                ax.pcolormesh(times, bazis, grad)

                ax.set_ylabel('back-azimuth')
                ax.set_xlabel('time')

                # xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
                # ax.xaxis.set_major_formatter(xfmt)
                # fig.autofmt_xdate()
                # fig.subplots_adjust(hspace=0)

                semblance = best_frame.reshape((n_bazis, n_slow))

                fig1 = self.pylab(name='FK: Max (%i)'%self.irun, get='figure')

                theta, r = num.meshgrid(bazis, slownesses)
                theta *= (num.pi/180.)

                ax = fig1.add_subplot(111, projection='polar')
                m = ax.pcolormesh(theta.T, r.T*km, to_db(semblance))

                self.adjust_polar_axis(ax)
                ax.plot(bazis[imax_bazi]*d2r, slownesses[imax_slow]*km, 'o')

                bazi_max = bazis[imax_bazi]*d2r
                slow_max = slownesses[imax_slow]*km
                ax.plot(bazi_max, slow_max, 'b.')
                ax.text(0.5, 0.01, 'Maximum at %s degrees, %s s/km' %
                        (num.round(bazi_max, 1), slow_max), transform=fig1.transFigure,
                        horizontalalignment='center', verticalalignment='bottom')

                if azi_theo:
                    ax.arrow(azi_theo/180.*num.pi, num.min(slownesses), 0,
                             num.max(slownesses), alpha=0.5, width=0.015,
                             edgecolor='black', facecolor='green', lw=2,
                             zorder=5)

                fig1.colorbar(m)

                fig1 = self.pylab(name='FK: Beam(%i)'%self.irun, get='figure')

                nsubplots = 4
                ax_raw = fig1.add_subplot(nsubplots, 1, 1)
                ax_shifted = fig1.add_subplot(nsubplots, 1, 2)
                ax_beam = fig1.add_subplot(nsubplots, 1, 3)

                axkwargs = dict(alpha=0.3, linewidth=0.3, color='grey')

                ybeam = num.zeros(lengthout)
                ybeam_weighted = num.zeros(lengthout)
                for i, (shift, array) in enumerate(zip(shifts.T, arrays)):
                    ax_raw.plot(times, array[npad: -npad], **axkwargs)
                    ishift = shift[imax_bazi_slow]
                    ax_shifted.plot(times, array[npad-ishift: -npad-ishift], **axkwargs)

                    ydata = traces[i].get_ydata()[npad-ishift: -npad-ishift]
                    ybeam += ydata

                    # calculate phase weighting
                    if self.want_pws:
                        ph_inst = instantaneous_phase(ydata)
                        ybeam_weighted += num.abs(num.exp(ph_inst))**4

                #ax_beam.plot(times, ybeam, color='black')
                ax_beam.plot(ybeam, color='black')
                ax_raw.set_title('Characteristic Function')
                ax_shifted.set_title('Shifted CF')
                ax_beam.set_title('Linear Stack')

                if self.want_pws:
                    ax_playground = fig1.add_subplot(nsubplots, 1, 4)
                    ax_playground.plot(ybeam*ybeam_weighted/len(arrays))
                    ax_playground.set_title('Phase Weighted Stack')

                beam_tr = trace.Trace(tmin=t_min+tpad, ydata=ybeam, deltat=deltat_cf)
                self.add_trace(beam_tr)

                self.irun += 1

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


def plot_back_slow_time(bazis, slownesses, times, frames, axs, window_length):
    '''
    :param window_length: in samples
    '''
    nbazis = len(bazis)
    nslow = len(slownesses)
    ntimes = len(times)

    #maxvals = num.argmax(frames, axis=0)
    maxvals = num.amax(frames, axis=0)
    # find maxima of *window_length* consecutive values:
    n = maxvals.size
    i_greatest_divisor = int(window_length/n)*n
    #imax = num.argmax(num.reshape(maxvals, (i_greatest_divisor, -1)), axis=1)

    # evaluate local extrema:
    #ilocal_max = argrelextrema(maxvals, num.greater)[0]

    #max_vs_t = num.argmax(frames[:, ilocal_max], axis=0)
    #ilocal_max_b, ilocal_max_s = num.unravel_index(max_vs_t, (nbazis, nslow))

    #ilocal_max_b, ilocal_max_s = num.unravel_index(imax, (nbazis, nslow))

    #max_times = times[ilocal_max]
    #max_semblance = maxvals[ilocal_max]
    #max_slownesses = slownesses[ilocal_max_s]
    #max_bazimuths = bazis[ilocal_max_b]
    #axs[0].plot(times, maxvals)
    axs[0].plot(maxvals)
    #axs[0].plot(max_times, max_semblance, 'x')
    #axs[1].scatter(max_times, max_slownesses, c=max_semblance, s=45)
    #axs[1].set_ylim((min(slownesses), max(slownesses)))
    #axs[2].scatter(max_times, max_bazimuths, c=max_semblance, s=45)


def __snufflings__():
    return [FK()]

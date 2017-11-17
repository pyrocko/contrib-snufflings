from __future__ import print_function
from builtins import range
import numpy as num
import logging, math

from pyrocko.gui.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko.gui.pile_viewer import Marker, EventMarker, PhaseMarker
from pyrocko.trace import Trace
from pyrocko import io, trace, util, cake, orthodrome, model
from pyrocko.dataset import crust2x2

km = 1000.
d2r = math.pi / 180.


class CorrelateEvents(Snuffling):

    def setup(self):
        self.set_name('Cross correlation relocation')
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 1.0,
            0.001, 50., low_is_none=True))
        self.add_parameter(Param('Lowpass [Hz]', 'corner_lowpass', 4.0,
            0.001, 50., high_is_none=True))
        self.add_parameter(Param('Time window begin', 'tstart', -1.0,
            -100., 0.))
        self.add_parameter(Param('Time window end', 'tend', 3.0,
            0., 100.))

        self.add_parameter(Param('Minimum correlation', 'min_corr', 0.5,
            0.0, 1.0))

        self.add_parameter(Param('Replace master depth [km]', 'master_depth_km', None,
            0.0, 100., high_is_none=True))

        self.add_parameter(Switch('Save figure', 'save', False))
        self.add_parameter(Switch('Fix depth', 'fix_depth', False))
        self.add_parameter(Switch('Show correlation traces', 'show_correlation_traces', False))
        self.add_parameter(Choice('Weighting', 'weighting', 'cubic',
            ['equal', 'linear', 'quadratic']))
        self.add_parameter(Choice('Earth model', 'model_select', 'Global',
            ['Global (ak135)', 'Local (from crust2x2)']))

        self.set_live_update(False)
        self.model = None
        self.model_key = None


    def call(self):
        self.cleanup()

        viewer = self.get_viewer()

        master = viewer.get_active_event()

        if master is None:
            self.fail('no master event selected')

        stations = list(viewer.stations.values())
        stations.sort(key=lambda s: (s.network,s.station))

        if not stations:
            self.fail('no station information available')

        # gather events to be processed

        events = []
        for m in viewer.markers:
            if isinstance(m, EventMarker):
                if m.kind == 0:
                    events.append( m.get_event() )

        events.sort(key=lambda ev: ev.time)

        event_to_number = {}
        for iev, ev in enumerate(events):
            event_to_number[ev] = iev

        if self.model_select.startswith('Global'):
            model_key = 'global'
        else:
            model_key = master.lat, master.lon

        if model_key != self.model_key:
            if self.model_select.startswith('Global'):
                self.model = cake.load_model()
            else:
                latlon = master.lat, master.lon
                profile = crust2x2.get_profile(*latlon)
                profile.set_layer_thickness(crust2x2.LWATER, 0.0)
                self.model = cake.LayeredModel.from_scanlines(
                        cake.from_crust2x2_profile(profile))

            self.model_key = model_key

        phases = {
                'P': ([ cake.PhaseDef(x) for x in 'P p'.split() ], 'Z'),
                'S': ([ cake.PhaseDef(x) for x in 'S s'.split() ], 'NE'),
            }

        phasenames = list(phases.keys())
        phasenames.sort()

        # synthetic arrivals and ray geometry for master event
        master_depth = master.depth
        if self.master_depth_km is not None:
            master_depth = self.master_depth_km * km

        tt = {}
        g = {}
        for iphase, phasename in enumerate(phasenames):
            for istation, station in enumerate(stations):
                dist = orthodrome.distance_accurate50m(master, station)
                azi = orthodrome.azimuth(master, station)

                arrivals = self.model.arrivals(
                        phases=phases[phasename][0],
                        distances=[ dist*cake.m2d ],
                        zstart = master_depth,
                        zstop = 0.0)

                if arrivals:
                    first = arrivals[0]
                    tt[station.network, station.station, phasename] = first.t

                    takeoff = first.takeoff_angle()
                    u = first.path.first_straight().u_in(first.endgaps)

                    g[iphase, istation] = num.array([
                            math.cos(azi*d2r) * math.sin(takeoff*d2r) * u,
                            math.sin(azi*d2r) * math.sin(takeoff*d2r) * u,
                            math.cos(takeoff*d2r) * u ])

        # gather picks for each event

        for ev in events:
            picks = {}
            for m2 in viewer.markers:
                if isinstance(m2, PhaseMarker) and m2.kind == 0:
                    if m2.get_event() == ev:
                        net, sta, _, _ = m2.one_nslc()
                        picks[net,sta,m2.get_phasename()] = (m2.tmax + m2.tmin) / 2.0

            ev.picks = picks

        # time corrections for extraction windows

        dataobs = []
        datasyn = []
        for phasename in phasenames:
            for station in stations:
                nsp = station.network, station.station, phasename
                datasyn.append(tt.get(nsp,None))
                for ev in events:
                    if nsp in ev.picks:
                        ttobs = ev.picks[nsp] - ev.time
                    else:
                        ttobs = None

                    dataobs.append(ttobs)

        ttsyn = num.array(datasyn, dtype=num.float).reshape((
            len(phasenames),
            len(stations)))

        ttobs = num.array(dataobs, dtype=num.float).reshape((
            len(phasenames),
            len(stations),
            len(events)))

        ttres = ttobs - ttsyn[:,:,num.newaxis]
        tt_corr_event = num.nansum( ttres, axis=1) /  \
                num.nansum( num.isfinite(ttres), axis=1 )

        tt_corr_event = num.where(num.isfinite(tt_corr_event), tt_corr_event, 0.)

        ttres -= tt_corr_event[:,num.newaxis,:]
        tt_corr_station = num.nansum( ttres, axis=2) /  \
                num.nansum( num.isfinite(ttres), axis=2 )

        tt_corr_station = num.where(num.isfinite(tt_corr_station), tt_corr_station, 0.)

        ttres -= tt_corr_station[:,:, num.newaxis]

        tevents_raw = num.array( [ ev.time for ev in events ] )

        tevents_corr = tevents_raw + num.mean(tt_corr_event, axis=0)

        # print timing information

        print('timing stats')

        for iphasename, phasename in enumerate(phasenames):
            data = []
            for ev in events:
                iev = event_to_number[ev]
                for istation, station in enumerate(stations):
                    nsp = station.network, station.station, phasename
                    if nsp in tt and nsp in ev.picks:
                        tarr = ev.time + tt[nsp]
                        tarr_ec = tarr + tt_corr_event[iphasename, iev]
                        tarr_ec_sc = tarr_ec + tt_corr_station[iphasename, istation]
                        tobs = ev.picks[nsp]

                        data.append((tobs-tarr, tobs-tarr_ec, tobs-tarr_ec_sc))

            if data:

                data = num.array(data, dtype=num.float).T

                print('event %10s %3s %3i %15.2g %15.2g %15.2g' % (
                        (ev.name, phasename, data.shape[1]) +
                            tuple( num.mean(num.abs(x)) for x in data )))
            else:
                print('event %10s %3s no picks' % (ev.name, phasename))

        # extract and preprocess waveforms

        tpad = 0.0
        for f in self.corner_highpass, self.corner_lowpass:
            if f is not None:
                tpad = max(tpad, 1.0/f)

        pile = self.get_pile()
        waveforms = {}
        for ev in events:
            iev = event_to_number[ev]
            markers = []
            for iphasename, phasename in enumerate(phasenames):
                for istation, station in enumerate(stations):
                    nsp = station.network, station.station, phasename
                    if nsp in tt:
                        tarr = ev.time + tt[nsp]
                        nslcs = [ ( station.network, station.station, '*', '*' ) ]
                        marker = PhaseMarker( nslcs, tarr, tarr, 1, event=ev,
                                phasename=phasename)
                        markers.append(marker)

                        tarr2 = tarr + tt_corr_station[iphasename, istation] + \
                                tt_corr_event[iphasename, iev]

                        marker = PhaseMarker( nslcs, tarr2, tarr2, 2, event=ev,
                                phasename=phasename)

                        markers.append(marker)

                        tmin = tarr2+self.tstart
                        tmax = tarr2+self.tend

                        marker = PhaseMarker( nslcs,
                                tmin, tmax, 3, event=ev,
                                phasename=phasename)

                        markers.append(marker)

                        trs = pile.all(tmin, tmax, tpad=tpad, trace_selector=
                                lambda tr: tr.nslc_id[:2] == nsp[:2],
                                want_incomplete=False)

                        trok = []
                        for tr in trs:
                            if num.all(tr.ydata[0] == tr.ydata):
                                continue

                            if self.corner_highpass:
                                tr.highpass(4, self.corner_highpass)
                            if self.corner_lowpass:
                                tr.lowpass(4, self.corner_lowpass)

                            tr.chop(tmin, tmax)
                            tr.set_location(ev.name)
                            #tr.shift( - (tmin - master.time) )

                            if num.all(num.isfinite(tr.ydata)):
                                trok.append(tr)

                        waveforms[nsp+(iev,)] = trok

            self.add_markers(markers)

        def get_channel(trs, cha):
            for tr in trs:
                if tr.channel == cha:
                    return tr
            return None

        nevents = len(events)
        nstations = len(stations)
        nphases = len(phasenames)

        # correlate waveforms

        coefs = num.zeros((nphases, nstations, nevents, nevents))
        coefs.fill(num.nan)
        tshifts = coefs.copy()
        tshifts_picked = coefs.copy()
        for iphase, phasename in enumerate(phasenames):
            for istation, station in enumerate(stations):
                nsp = station.network, station.station, phasename

                for a in events:
                    ia = event_to_number[a]
                    for b in events:
                        ib = event_to_number[b]

                        if ia == ib:
                            continue

                        if nsp in a.picks and nsp in b.picks:
                            tshifts_picked[iphase,istation,ia,ib] = \
                                    b.picks[nsp] - a.picks[nsp]

                        wa = waveforms[nsp+(ia,)]
                        wb = waveforms[nsp+(ib,)]

                        channels = list(set([ tr.channel for tr in wa + wb ]))
                        channels.sort()

                        tccs = []
                        for cha in channels:
                            if cha[-1] not in phases[phasename][1]:
                                continue

                            ta = get_channel(wa, cha)
                            tb = get_channel(wb, cha)
                            if ta is None or tb is None:
                                continue

                            tcc = trace.correlate(ta,tb, mode='full', normalization='normal',
                                    use_fft=True)

                            tccs.append(tcc)

                        if not tccs:
                            continue

                        tc = None
                        for tcc in tccs:
                            if tc is None:
                                tc = tcc
                            else:
                                tc.add(tcc)

                        tc.ydata *= 1./len(tccs)

                        tmid = tc.tmin*0.5 + tc.tmax*0.5
                        tlen = (tc.tmax - tc.tmin)*0.5
                        tc_cut = tc.chop(tmid-tlen*0.5, tmid+tlen*0.5, inplace=False)

                        tshift, coef = tc_cut.max()

                        if (tshift < tc.tmin + 0.5*tc.deltat or
                                tc.tmax - 0.5*tc.deltat < tshift):
                            continue

                        coefs[iphase,istation,ia,ib] = coef
                        tshifts[iphase,istation,ia,ib] = tshift

                        if self.show_correlation_traces:
                            tc.shift(master.time - (tc.tmax + tc.tmin)/2.)
                            self.add_trace(tc)


        #tshifts = tshifts_picked

        coefssum_sta = num.nansum(coefs, axis=2) / num.sum(num.isfinite(coefs), axis=2)
        csum_sta = num.nansum(coefssum_sta, axis=2) / num.sum(num.isfinite(coefssum_sta), axis=2)

        for iphase, phasename in enumerate(phasenames):
            for istation, station in enumerate(stations):
                print('station %-5s %s %15.2g' %
                      (station.station, phasename, csum_sta[iphase,istation]))

        coefssum = num.nansum(coefs, axis=1) / num.sum(num.isfinite(coefs), axis=1)
        csumevent = num.nansum(coefssum, axis=2) / num.sum(num.isfinite(coefssum), axis=2)
        above = num.where(num.isfinite(coefs), coefs >= self.min_corr, 0)

        csumabove = num.sum(num.sum(above, axis=1), axis=2)

        coefssum = num.ma.masked_invalid(coefssum)

        print('correlation stats')

        for iphase, phasename in enumerate(phasenames):
            for ievent, event in enumerate(events):
                print('event %10s %3s %8i %15.2g' % (
                        event.name, phasename,
                        csumabove[iphase,ievent], csumevent[iphase,ievent]))

        # plot event correlation matrix

        fframe = self.figure_frame()
        fig = fframe.gcf()

        for iphase, phasename in enumerate(phasenames):

            p = fig.add_subplot(1,nphases,iphase+1)
            p.set_xlabel('Event number')
            p.set_ylabel('Event number')
            mesh = p.pcolormesh(coefssum[iphase])
            cb = fig.colorbar(mesh, ax=p)
            cb.set_label('Max correlation coefficient')

        if self.save:
            fig.savefig(self.output_filename(dir='correlation.pdf'))

        fig.canvas.draw()


        # setup and solve linear system

        data = []
        rows = []
        weights = []
        for iphase in range(nphases):
            for istation in range(nstations):
                for ia in range(nevents):
                    for ib in range(ia+1,nevents):
                        k = iphase, istation, ia, ib
                        w = coefs[k]
                        if not num.isfinite(tshifts[k]) \
                                or not num.isfinite(w) or w < self.min_corr:
                            continue

                        row = num.zeros(nevents*4)
                        row[ia*4:ia*4+3] = g[iphase,istation]
                        row[ia*4+3] = -1.0
                        row[ib*4:ib*4+3] = -g[iphase,istation]
                        row[ib*4+3] = 1.0
                        weights.append(w)

                        rows.append(row)
                        data.append(tshifts[iphase,istation,ia,ib])

        nsamp = len(data)

        for i in range(4):
            row = num.zeros(nevents*4)
            row[i::4] = 1.
            rows.append(row)
            data.append(0.0)

        if self.fix_depth:
            for ievent in range(nevents):
                row = num.zeros(nevents*4)
                row[ievent*4+2] = 1.0
                rows.append(row)
                data.append(0.0)

        a = num.array(rows, dtype=num.float)
        d = num.array(data, dtype=num.float)
        w = num.array(weights, dtype=num.float)

        if self.weighting == 'equal':
            w[:nsamp] = 1.0
        elif self.weighting == 'linear':
            pass
        elif self.weighting == 'quadratic':
            w[:nsamp] = w[:nsamp]**2

        a[:nsamp,:] *= w[:,num.newaxis]
        d[:nsamp] *= w[:nsamp]

        x, residuals, rank, singular = num.linalg.lstsq(a,d)

        x0 = num.zeros(nevents*4)
        x0[3::4] = tevents_corr
        mean_abs_residual0 = num.mean(
                num.abs((num.dot(a[:nsamp], x0) - d[:nsamp])/w[:nsamp]))

        mean_abs_residual = num.mean(
                num.abs((num.dot(a[:nsamp],x) - d[:nsamp])/w[:nsamp]))

        print(mean_abs_residual0, mean_abs_residual)

        # distorted solutions

        npermutations = 100
        noiseamount = mean_abs_residual
        xdistorteds = []
        for i in range(npermutations):
            dnoisy = d.copy()
            dnoisy[:nsamp] += num.random.normal(size=nsamp)*noiseamount*w[:nsamp]
            xdistorted, residuals, rank, singular = num.linalg.lstsq(a,dnoisy)
            xdistorteds.append(xdistorted)

            mean_abs_residual = num.mean(num.abs(num.dot(a,xdistorted)[:nsamp] - dnoisy[:nsamp]))

        tmean = num.mean([ e.time for e in events ])

        north = x[0::4]
        east = x[1::4]
        down = x[2::4]
        etime = x[3::4] + tmean

        def plot_range(x):
            mi, ma = num.percentile(x, [10., 90.])
            ext = (ma-mi)/5.
            mi -= ext
            ma += ext
            return mi, ma

        lat, lon = orthodrome.ne_to_latlon(master.lat, master.lon, north, east)

        events_out = []
        for ievent, event in enumerate(events):
            event_out = model.Event(time=etime[ievent],
                    lat=lat[ievent],
                    lon=lon[ievent],
                    depth=down[ievent] + master_depth,
                    name = event.name)

            mark = EventMarker(event_out, kind=4)
            self.add_marker(mark)
            events_out.append(event_out)

        model.Event.dump_catalog(events_out, 'events.relocated.txt')

        # plot results

        ned_orig = []
        for event in events:
            n, e = orthodrome.latlon_to_ne(master, event)
            d = event.depth

            ned_orig.append((n,e,d))

        ned_orig = num.array(ned_orig)

        ned_orig[:,0] -= num.mean(ned_orig[:,0])
        ned_orig[:,1] -= num.mean(ned_orig[:,1])
        ned_orig[:,2] -= num.mean(ned_orig[:,2])

        north0, east0, down0 = ned_orig.T

        north2, east2, down2, time2 = num.hstack(xdistorteds).reshape((-1,4)).T

        fframe = self.figure_frame()
        fig = fframe.gcf()

        color_sym = (0.1,0.1,0.0)
        color_scat = (0.3,0.5,1.0,0.2)

        d = u'\u0394 '

        if not self.fix_depth:
            p = fig.add_subplot(2,2,1, aspect=1.0)
        else:
            p = fig.add_subplot(1,1,1, aspect=1.0)

        mi_north, ma_north = plot_range(north)
        mi_east, ma_east = plot_range(east)
        mi_down, ma_down = plot_range(down)

        p.set_xlabel(d+'East [km]')
        p.set_ylabel(d+'North [km]')
        p.plot(east2/km, north2/km, '.', color=color_scat, markersize=2)
        p.plot(east/km, north/km, '+', color=color_sym)
        p.plot(east0/km, north0/km, 'x', color=color_sym)
        p0 = p

        for i,ev in enumerate(events):
            p.text(east[i]/km, north[i]/km, ev.name, clip_on=True)

        if not self.fix_depth:

            p = fig.add_subplot(2,2,2, sharey=p0, aspect=1.0)
            p.set_xlabel(d+'Depth [km]')
            p.set_ylabel(d+'North [km]')
            p.plot(down2/km, north2/km, '.', color=color_scat, markersize=2)
            p.plot(down/km, north/km, '+', color=color_sym)
            for i,ev in enumerate(events):
                p.text(down[i]/km, north[i]/km, ev.name, clip_on=True)

            p1 = p

            p = fig.add_subplot(2,2,3, sharex=p0, aspect=1.0)
            p.set_xlabel(d+'East [km]')
            p.set_ylabel(d+'Depth [km]')
            p.plot(east2/km, down2/km, '.', color=color_scat, markersize=2)
            p.plot(east/km, down/km, '+', color=color_sym)
            for i,ev in enumerate(events):
                p.text(east[i]/km, down[i]/km, ev.name, clip_on=True)


            p.invert_yaxis()
            p2 = p

        p0.set_xlim(mi_east/km, ma_east/km)
        p0.set_ylim(mi_north/km, ma_north/km)

        if not self.fix_depth:
            p1.set_xlim(mi_down/km, ma_down/km)
            p2.set_ylim(mi_down/km, ma_down/km)

        if self.save:
            fig.savefig(self.output_filename(dir='locations.pdf'))

        fig.canvas.draw()


def __snufflings__():
    return [CorrelateEvents()]


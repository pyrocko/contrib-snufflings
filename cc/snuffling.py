import sys
import os.path as op
import numpy as num
import progressbar

from pyrocko import model, trace, pile, util, orthodrome, cake, gui_util
from pyrocko.gf.seismosizer import Target, SeismosizerTrace
from pyrocko.snuffling import Snuffling, Choice, Param, Switch
from similarity import SimilarityMatrix, Similarity

util.setup_logging('cc.py')


def make_targets(pile, stations, channels):
    targets = []
    for nslc_id in pile.nslc_ids.keys():
        for s in stations:
            if util.match_nslc('%s.*%s'%(s.nsl_string(), channels), nslc_id):
                targets.append(Target(lat=s.lat,
                                     lon=s.lon,
                                     depth=s.depth,
                                     elevation=s.elevation,
                                     codes=nslc_id))
            else:
                continue
    return targets


class SimilaritySnuffling(Snuffling):
    '''
    tdist: Time distance from center of cross correlated trace. Sets the intended time span.
    vmin and vmax: horizontal velocities used to chop the traces to be cross correlated
    low and high: low- and high pass corner frequencies of 4th order butterworth filter
    '''
    def setup(self):
        self.set_name('Similarity')
        self.add_parameter(Choice('Channels', 'channels', 'Z', ['Z', '*']))
        self.add_parameter(Param('limit', 'limit', 0.9, 0.0, 1.0))
        # Todo visibility changed -> hp lp
        self.add_parameter(Param('lp', 'low', 10., 0.1, 200.0, high_is_none=True))
        self.add_parameter(Param('hp', 'high', 1, 0.1, 200.0, low_is_none=True))
        self.add_parameter(Param('padding', 'tpad', 10, 0.1, 60.0))
        self.add_parameter(Param('dt wanted', 'dt_wanted', 0.01, 0.01, 10., low_is_none=True))
        self.add_parameter(Param('tdist', 'tdist', 7.5,1., 20.))
        self.add_parameter(Param('v min', 'vmin', 1500., 500., 2000.))
        self.add_parameter(Param('v max', 'vmax', 6000., 2000., 1000.))
        self.add_parameter(Switch('show traces', 'show_traces', False))
        self.add_parameter(Switch('show results', 'show_results', False))
        self.add_parameter(Switch('Apply to full dataset', 'apply_full', False))
        self.add_trigger('save result', self.save)
        self.set_live_update(False)

    def call(self):
        self.cleanup()
        events = [m.get_event() for m in self.get_selected_event_markers()]
        for iev, ev in enumerate(events):
            ev.name = '%05i' % iev


        show_arrivals = False

        filters = []
        for ident in ['high', 'low']:
            val = getattr(self, ident)
            if val != None:
                filters.append(trace.ButterworthResponse(corner=float(val),
                                                         order=4,
                                                         type=ident))

        stations = self.get_stations()
        stations_d = dict((s.nsl(), s) for s in stations)

        # TODO option to choose other models
        mod = cake.load_model()
        nevents = len(events)
        pile = self.get_pile()
        targets = make_targets(pile, stations, self.channels)
        if len(targets)==0:
            self.fail("No station available")
        ntargets = len(targets)
        self.cc = num.zeros((ntargets, nevents, nevents), dtype=num.float)
        self.similarity_matrix = SimilarityMatrix(targets=targets,
                                                  events=events,
                                                  filters=filters,
                                                  padding=float(self.tpad))
        similarities = []
        trs2add = []
        if self.apply_full:
            tmin_selection = None; tmax_selection=None
        else:
            tmin_selection, tmax_selection = self.get_selected_time_range(fallback=True)
        istep = 0
        #pb = progressbar.ProgressBar(maxval=ntargets).start()
        for itarget, target in enumerate(targets):
            ok_raw = []
            ok_filtered = []
            markers = []
            for iev, ev in enumerate(events):
                dist = target.distance_to(ev)
                tmin = ev.time + dist / self.vmax - self.tpad
                tmax = ev.time + dist / self.vmin + self.tpad
                trs = pile.chopper(tmin=tmin,
                                   tmax=tmax,
                                   trace_selector=lambda tr: tr.nslc_id == target.codes,
                                   want_incomplete=False)
                tr = [t for trss  in trs for t in trss]

                if len(tr)==0:
                    continue
                elif len(tr)==1:
                    tr = tr[0]
                else:
                    self.fail('Something went wrong')
                if self.dt_wanted:
                    tr.downsample_to(self.dt_wanted)

                tr2 = tr.copy()
                for f in filters:
                    tr2.transfer(transfer_function=f)
                tr2.chop(tmin, tmax)
                tr2.set_codes(location=ev.name+'f')

                tr.chop(tmin, tmax)
                tr.set_codes(location=ev.name+'r')

                #ok_raw.append((iev, ev, tr))
                ok_filtered.append((iev, ev, tr2))

                #if show_arrivals:
                #    rays = mod.arrivals(
                #        phases=[cake.PhaseDef(x) for x in 'p P'.split()],
                #        distances=[dist*cake.m2d], 
                #        zstart=ev.depth)

                #    if rays:
                #        ray = rays[0]
                #        mark = gui_util.PhaseMarker(
                #            [tr.nslc_id],
                #            ev.time + ray.t,
                #            ev.time + ray.t,
                #            0,
                #            phasename='P')

                #        #print mark

                #        markers.append(mark)


            ok = ok_filtered
            # isn't it doing the same thing twice?
            # better to pop items?
            while ok:
                (ia, a_ev, a_tr) = ok.pop()
                #for (ia, a_ev, a_tr) in ok:
                for (ib, b_ev, b_tr) in ok:
                    relamp = 0.0
                    if a_tr is not None and b_tr is not None:
                        c_tr = trace.correlate(a_tr, b_tr, mode='full', normalization='normal')
                        t_center = c_tr.tmin+(c_tr.tmax-c_tr.tmin)/2.
                        c_tr_chopped = c_tr.chop(t_center-self.tdist, t_center+self.tdist, inplace=False)
                        t_mini, v_mini = c_tr_chopped.min()
                        t_maxi, v_maxi = c_tr_chopped.max()
                        #print t_mini, t_maxi, v_mini, v_maxi, (abs(v_mini) > abs(v_maxi) and v_mini < 0.95)

                        b_tr_shifted = b_tr.copy()
                        a_tr_shifted = a_tr.copy()

                        # TODO: Confirm that traces are shifted correctly!!!!!
                        if abs(v_mini) > abs(v_maxi):
                            if abs(t_center-t_mini) < self.limit * self.tdist:
                                self.cc[itarget, ia, ib] = v_mini
                                b_tr_shifted.shift(-t_mini)
                        else:
                            if abs(t_center-t_maxi) < self.limit * self.tdist:
                                self.cc[itarget, ia, ib] = v_maxi
                                b_tr_shifted.shift(-t_maxi)

                        if self.cc[itarget, ia, ib] != 0.0:
                            tmin = max(a_tr.tmin, b_tr_shifted.tmin)
                            tmax = min(a_tr.tmax, b_tr_shifted.tmax)
                            a_tr_chopped = a_tr.chop(tmin, tmax, inplace=False)
                            b_tr_chopped = b_tr_shifted.chop(tmin, tmax)

                            ya = a_tr_chopped.ydata
                            yb = b_tr_chopped.ydata
                            relamp = num.sum(ya*yb) / num.sum(ya**2)

                        if self.show_traces:
                            #trs2add.extend([tr.copy() for (iev, ev, tr) in (ok_raw + ok_filtered) if tr])
                            trs2add.extend([c_tr])
                            #for tr in trs2add:
                            #tr.set_location('%s%s'%(tr.location, 'cc'))
                        sim = Similarity(
                            ievent=ia,
                            jevent=ib,
                            itarget=itarget,
                            cross_correlation=float(self.cc[itarget, ia, ib]),
                            cross_correlation_trace=SeismosizerTrace.from_pyrocko_trace(c_tr),
                            relative_amplitude=float(relamp))

                        similarities.append(sim)

        if self.show_results:
            for itarget, target in enumerate(targets):
                fig = self.pylab(get='figure')
                axes = fig.add_subplot(111)
                axes.set_xlabel('Event number')
                axes.set_ylabel('Event number')
                mesh = axes.pcolormesh(self.cc[itarget,:,:], cmap='RdBu', vmin=-1.0, vmax=1.0)
                cb = fig.colorbar(mesh, ax=axes)
                cb.set_label('Max correlation coefficient')
                fig.canvas.draw()

        if self.show_traces:
            self.add_traces(trs2add)
        self.similarity_matrix.similarities = similarities
        self.similarity_matrix.validate()

    def save(self):
        output_filename = self.output_filename()
        self.similarity_matrix.dump(filename=output_filename)

def __snufflings__():
    return [SimilaritySnuffling()]

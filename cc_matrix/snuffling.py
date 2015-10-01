import sys
import matplotlib
matplotlib.use('Qt4Agg')
import logging
import os.path as op
import numpy as num
import progressbar

from pyrocko import model, trace, util, orthodrome, cake, gui_util
from pyrocko.gf.seismosizer import Target, SeismosizerTrace
from pyrocko.snuffling import Snuffling, Choice, Param, Switch
from similarity import SimilarityMatrix, Similarity
import matplotlib.pyplot as plt
util.setup_logging('cc.py')
logger = logging.getLogger('cc-snuffling')

def make_targets(pile, stations):
    targets = []
    for nslc_id in pile.nslc_ids.keys():
        for s in stations:
            if util.match_nslc('%s.*'%(s.nsl_string()), nslc_id):
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
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Cross Correlation Matrix</h1>
    <p>
    Cross correlate events and save results in yaml format. Requires station
    meta information to be availble. Selected events will be cross correlated
    using data from visible traces, exclusively.
    <b>Parameters:</b><br />
        <b>&middot; Windowing method:</b> Time window selection behaviour.
        Select <i>P-phase</i> to use P phase onset to set the position of the
        selection window.<br>
        if windowing method is using P-phase, the chopping ranges are caluclated as follows:
        <p>
        tmin = t_p-padding*0.1
        <br>
        tmax = t_p+padding*0.9
        <p>
        where t_p is the earlier onset time of the p- and P-phase respectively.<br>
        Select <i>vmin/vmax</i> and the use sliders <i>vmin</i> and <i>vmax</i> to
        define horizontal velocities to define the selection window.<br />
        <b>&middot; tdist:</b> Time distance from center of cross correlated
        trace. Sets the intended time span. <br />
        <b>&middot; dt wanted:</b> downsample traces to selected delta t.<br />
        <b>&middot; low and high:</b> Low- and high pass corner frequencies of
        4th order butterworth filter.<br />
        <b>&middot; save traces:</b> Select an output directory using the opening
        dialog. Trace pairs, as well as the cc trace are going to be saved here.
        <b>Drawing the traces is time consuming!</b><br />
        <b>&middot; show results:</b> make cc result images after processing.<br />
        <b>&middot; save results:</b> Store results in YAML format.<br />
        <p>
        Since results are stored in yaml format, they can easily be loaded
        programmatically as follows:
        <pre>
        from similarity import SimilarityMatrix, Similarity
        from pyrocko.guts import load as guts_load

        fn = 'saved-similarities.dat'
        matrix = guts_load(filename=fn)
        for s in matrix.similarities:
                print s
        </pre>
        This assumes, that you the module <pre>similarity.py</pre> is stored in
        the same directory or in a place where it can be found by python
        e.g. using the the <pre>PYTHONPATH</pre> environment variable.
    </p>
    '''
    def setup(self):
        self.set_name('CC Matrix')
        self.add_parameter(Choice(
            'Windowing method', 'time_window_choice', 'P-phase',
            ['P-phase', 'vmin/vmax']))
        self.add_parameter(Param('low', 'low', 10., 0.1, 200.0,
                                 high_is_none=True))
        self.add_parameter(Param('high', 'high', 1, 0.1, 200.0,
                                 low_is_none=True))
        self.add_parameter(Param('padding [s]', 'tpad', 10, 0.1, 60.0))
        self.add_parameter(Param('dt wanted', 'dt_wanted', 0.01, 0.01, 10.,
                                 low_is_none=True))
        self.add_parameter(Param('tdist [s]', 'tdist', 7.5,1., 20.))
        self.add_parameter(Param('v min [m/s]', 'vmin', 1500., 500., 2000.))
        self.add_parameter(Param('v max [m/s] ', 'vmax', 2000., 6000., 1000.))
        self.add_parameter(Switch('Save Traces', 'save_traces', False))
        self.add_parameter(Switch('Show Results', 'show_results', False))
        self.add_trigger('Save Result', self.save)
        self.set_live_update(False)
        self.phase_cache = {}

    def call(self):
        self.cleanup()
        viewer = self.get_viewer()
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
        traces = list(self.chopper_selected_traces(fallback=True, trace_selector=
                                                   viewer.trace_selector,
                                                   load_data=False))
        traces = [tr for trs in traces for tr in trs ]
        visible_nslcs = [tr.nslc_id for tr in traces]
        stations = filter(lambda s: util.match_nslcs("%s.%s.%s.*" % s.nsl(), visible_nslcs), stations)

        # TODO option to choose other models
        mod = cake.load_model()
        nevents = len(events)
        pile = self.get_pile()
        targets = make_targets(pile, stations)
        if len(targets)==0:
            self.fail("No station available")
        ntargets = len(targets)
        self.cc = num.zeros((ntargets, nevents, nevents), dtype=num.float)
        self.similarity_matrix = SimilarityMatrix(targets=targets,
                                                  events=events,
                                                  filters=filters,
                                                  padding=float(self.tpad),
                                                  windowing_method=self.time_window_choice,
                                                  vmax=float(self.vmax),
                                                  vmin=float(self.vmin))
        similarities = []
        trs2add = []
        if self.save_traces :
            figure_dir = self.input_directory(caption='Select directory to store images')
        for itarget, target in enumerate(targets):
            print (itarget+1.)/float(ntargets)
            ok_filtered = []
            markers = []
            for iev, ev in enumerate(events):
                dist = target.distance_to(ev)
                if self.time_window_choice=='vmin/vmax':
                    tmin = ev.time + dist / self.vmax - self.tpad
                    tmax = ev.time + dist / self.vmin + self.tpad
                elif self.time_window_choice=='P-phase':
                    d = dist*cake.m2d
                    z = ev.depth
                    t = self.phase_cache.get((mod, d, z), False)
                    if not t:
                        rays = mod.arrivals(
                            phases=[cake.PhaseDef(x) for x in 'p P'.split()],
                            distances=[d],
                            zstart=z)
                        t = rays[0].t
                        self.phase_cache[(mod, d, z)] = t
                    tmin = ev.time + t - self.tpad * 0.1
                    tmax = ev.time + t + self.tpad * 0.9
                trs = pile.chopper(tmin=tmin,
                                   tmax=tmax,
                                   trace_selector=viewer.trace_selector,
                                   want_incomplete=False)
                tr = [t for trss  in trs for t in trss if t.nslc_id==target.codes]
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
                    tr2 = tr2.transfer(transfer_function=f)
                tr2.chop(tmin, tmax)
                tr2.set_codes(location=ev.name+'f')

                tr.chop(tmin, tmax)
                tr.set_codes(location=ev.name+'r')

                ok_filtered.append((iev, ev, tr2))

            ok = ok_filtered
            while ok:
                (ia, a_ev, a_tr) = ok.pop()
                for (ib, b_ev, b_tr) in ok:
                    relamp = 0.0
                    if a_tr is not None and b_tr is not None:
                        c_tr = trace.correlate(a_tr, b_tr, mode='full', normalization='normal')
                        t_center = c_tr.tmin+(c_tr.tmax-c_tr.tmin)/2.
                        c_tr_chopped = c_tr.chop(t_center-self.tdist, t_center+self.tdist, inplace=False)
                        t_mini, v_mini = c_tr_chopped.min()
                        t_maxi, v_maxi = c_tr_chopped.max()
                        b_tr_shifted = b_tr.copy()
                        a_tr_shifted = a_tr.copy()

                        if abs(v_mini) > abs(v_maxi):
                            v_cc = v_mini
                            time_lag = -t_mini
                        else:
                            time_lag = -t_maxi
                            v_cc = v_maxi

                        self.cc[itarget, ia, ib] = v_cc
                        b_tr_shifted.shift(time_lag)

                        if self.cc[itarget, ia, ib] != 0.0:
                            tmin = max(a_tr.tmin, b_tr_shifted.tmin)
                            tmax = min(a_tr.tmax, b_tr_shifted.tmax)
                            try:
                                a_tr_chopped = a_tr.chop(tmin, tmax, inplace=False)
                                b_tr_chopped = b_tr_shifted.chop(tmin, tmax)
                            except trace.NoData:
                                logger.warn('NoData %s'%a_tr_chopped)
                                continue

                            ya = a_tr_chopped.ydata
                            yb = b_tr_chopped.ydata
                            relamp = num.sum(ya*yb) / num.sum(ya**2)

                        if self.save_traces:
                            fig, axes = plt.subplots(3,1)
                            fig.suptitle('.'.join(target.codes))
                            axes[0].plot(a_tr_chopped.get_xdata(), a_tr_chopped.get_ydata())
                            axes[0].text(0, 1, "id: %s, time: %s" %(a_ev.name, util.time_to_str(a_ev.time)),
                                         transform=axes[0].transAxes,
                                         verticalalignment='top', horizontalalignment='left')
                            axes[1].plot(b_tr_chopped.get_xdata(), b_tr_chopped.get_ydata())
                            axes[1].text(0, 1, "id: %s, time: %s" %(b_ev.name, util.time_to_str(b_ev.time)),
                                         transform=axes[1].transAxes,
                                         verticalalignment='top', horizontalalignment='left')
                            axes[2].plot(c_tr.get_xdata(), c_tr.get_ydata())
                            axes[2].text(0, 1, 'cc_max: %1.4f' % v_cc,
                                         transform=axes[2].transAxes,
                                         verticalalignment='top', horizontalalignment='left')
                            fn = op.join(figure_dir, 'cc_T%s.E%s.E%s.png' % (itarget, ia, ib))
                            fig.savefig(fn, pad_inches=0.1, bbox_inches='tight', tight_layout=True)

                        sim = Similarity(
                            ievent=ia,
                            jevent=ib,
                            itarget=itarget,
                            cross_correlation=float(self.cc[itarget, ia, ib]),
                            relative_amplitude=float(relamp),
                            time_lag=float(-time_lag))

                        similarities.append(sim)

        if self.show_results:
            for itarget, target in enumerate(targets):
                if not num.any(self.cc[itarget]):
                    continue
                fig = self.pylab(get='figure')
                fig.suptitle('.'.join(target.codes))
                axes = fig.add_subplot(111)
                axes.set_xlabel('Event number')
                axes.set_ylabel('Event number')
                mesh = axes.pcolormesh(self.cc[itarget,:,:], cmap='RdBu', vmin=-1.0, vmax=1.0)
                cb = fig.colorbar(mesh, ax=axes)
                cb.set_label('Max correlation coefficient')
                fig.canvas.draw()

        self.similarity_matrix.similarities = similarities
        self.similarity_matrix.validate()

    def save(self):
        output_filename = self.output_filename()
        self.similarity_matrix.dump(filename=output_filename)

def __snufflings__():
    return [SimilaritySnuffling()]

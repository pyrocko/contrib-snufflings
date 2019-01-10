from pyrocko.gui.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko.gui.pile_viewer import Marker, EventMarker, PhaseMarker
from pyrocko import io, trace, model
from collections import defaultdict


class CorrsearchSnuffling(Snuffling):

    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h2 align="center">Cross Correlation Pattern Matching</h2>
    <h3 align="left">Usage</h3>
    <p>Select master signal using extended markers. Only one master signal per
    channel can be used.</p>
    <p>
    <b>Set parameters:</b><br />
        <b>&middot; Downsample to [Hz]</b>  -  Reduce the number of samples to
                improve speed<br />
        <b>&middot; Highpass[Hz]/Lowpass[Hz]</b>  -  Filter traces before cross
                correlation<br />
        <b>&middot; Apply to full dataset</b>  -  Work on the entire data set
                (By default, processing is limited to data shown in viewer)
                <br />
        <b>&middot; Normalization</b>  -  Normalize traces before cross
                correlation. <i>Gliding</i>: time varying normalization (using
                moving sum); <i>Normal</i>: normalizes by product of energy of
                both traces<br />
        <b>&middot; Threshold</b>  -  Peak detection threshold applied to
                correlation trace<br />
    </p>
    Hit <b>run</b>.<br>
    Cross correlations of different channels are summed and added to the viewer
    as an additional trace labeled <i>Sum Cross Correlation</i>.
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Cross Correlation Search')
        self.add_parameter(
            Param('Downsample to [Hz]', 'downsample',
            None, 0.1, 200., high_is_none=True))
        
        self.add_parameter(
            Param('Highpass [Hz]', 'corner_highpass',
            None, 0.001, 50., low_is_none=True))
        
        self.add_parameter(
            Param('Lowpass [Hz]', 'corner_lowpass',
            None, 0.001, 50., high_is_none=True))
        
        self.add_parameter(Param('t search', 'tsearch', 1., 0.1, 10.))
        self.add_parameter(Switch('Apply to full dataset', 'apply_to_all', False))
        self.add_parameter(Choice('Normalization', 'normalization', 'Off',
            ('Off', 'Normal', 'Gliding')))
        self.add_parameter(Param('Threshold', 'threshold', 0.5, 0., 1.))
        self.add_parameter(Switch('Use FFT', 'use_fft', True))
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()

        if self.corner_highpass:
            tpad = 1./self.corner_highpass
        else:
            tpad = 0.

        try:
            viewer = self.get_viewer()
            markers = viewer.selected_markers()
            if not markers:
                return

            if len(markers) != 1:
                return

            marker = markers[0]
            master_tmin, master_tmax = marker.tmin, marker.tmax
            if master_tmin >= master_tmax:
                return

        except NoViewerSet:
            viewer = None
            master_tmin, master_tmax = self.master_tmin, self.master_tmax

        pile = self.get_pile()
        masters = {}

        def preprocess(_tr):
            if self.downsample:
                _tr.downsample_to(1./self.downsample)

            if self.corner_highpass:
                _tr.highpass(4, self.corner_highpass)

            if self.corner_lowpass:
                _tr.lowpass(4, self.corner_lowpass)

        for tr in pile.all(tmin=master_tmin, tmax=master_tmax, tpad=tpad):
            for m in markers:
                if m.match_nslc(tr.nslc_id):
                    preprocess(tr)
                    tr.chop(tr.wmin, tr.wmax)
                    masters[tr.nslc_id] = tr
                    break

        if self.apply_to_all:
            tmin, tmax = pile.get_tmin()+tpad, pile.get_tmax()
        else:
            tmin, tmax = self.get_viewer().get_time_range()

        normalization = {'Off': None, 'Normal': 'normal', 'Gliding': 'gliding'}[self.normalization]

        for traces in pile.chopper(tmin=tmin, tmax=tmax, want_incomplete=True):
            sccs = defaultdict()
            sccn = defaultdict()
            for b in traces:
                nslc = b.nslc_id
                if nslc in masters:
                    a = masters[nslc]
                    preprocess(b)

                    c = trace.correlate(
                        a, b, mode='valid',
                        normalization=normalization,
                        use_fft=self.use_fft)

                    c.shift(-c.tmin + b.tmin)
                    c.meta = {'tabu': True}

                    scc = sccs.get(nslc, None)
                    if not scc:
                        scc = c.copy()
                        scc.meta = {'tabu': True}
                        scc.wmin = b.wmin
                        scc.wmax = b.wmax
                        scc.set_codes(location=scc.location+'_SUM')
                        sccn[nslc] = 1
                    else:
                        sccn[nslc] += 1

                    sccs[nslc] = scc

            for nslc_id, scc in sccs.items():
                scc.ydata /= sccn[nslc_id]
                scc.chop(scc.wmin, scc.wmax)

                markers = []
                for t, a in zip(*scc.peaks(self.threshold, tsearch=self.tsearch)):
                    m = PhaseMarker(tmin=t, tmax=t, phasename='%1.3f' % a, kind=3, nslc_ids=(nslc_id,))
                    markers.append(m)

                if viewer:
                    self.add_traces([scc])
                    self.add_markers(markers)
                else:
                    io.save([scc], self.out_path, format='from_extension')


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''

    return [CorrsearchSnuffling()]


if __name__ == '__main__':

    snuf = CorrsearchSnuffling()
    snuf.setup()
    snuf.apply_to_all = True
    snuf.corner_highpass = 0.1
    markers = Marker.load_markers('event.picks')
    m = markers[0]
    snuf.master_tmin, snuf.master_tmax = m.tmin, m.tmax
    snuf.out_path = 'corr/%(tmin)s.yaff'
    snuf.call()

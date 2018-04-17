from pyrocko.gui.snuffling import Snuffling, Param, Switch
import numpy as num

from pyrocko.plot import graph_colors as colors


def to01(c):
    return c[0]/255., c[1]/255., c[2]/255.


class PlotPSD(Snuffling):

    '''
    <html>
    <body>
    <h1>Plot PSD (Power Spectral Density)</h1>

    Visible or selected data is cut into windows of 2 x 'Window length',
    tapered with a Hanning taper, FFTed, sqared, normalized and gathered by
    mean or median and percentiles.

    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Plot PSD')
        self.add_parameter(
            Param('Window length [s]:', 'tinc', 100, 0.1, 10000.,
                  high_is_none=True))
        self.add_parameter(Switch('Save figure', 'save', False))
        self.add_parameter(Switch('Join stations', 'join_stations', False))
        self.add_parameter(Switch('Show mean', 'mean', False))
        self.add_parameter(Switch('Show logmean', 'logmean', False))
        self.add_parameter(Switch('Show median', 'median', True))
        self.add_parameter(Switch('Show percentiles', 'percentiles', False))
        self.add_parameter(Switch('Show min and max', 'minmax', False))
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''

        by_nslc = {}
        if self.tinc is not None:
            tpad = self.tinc/2
        else:
            tpad = 0.0

        for traces in self.chopper_selected_traces(
                tinc=self.tinc, tpad=tpad,
                want_incomplete=False, fallback=True):
            for tr in traces:
                nslc = tr.nslc_id
                if self.tinc is not None:
                    nwant = int(self.tinc * 2 / tr.deltat)

                    if nwant != tr.data_len():
                        if tr.data_len() == nwant + 1:
                            tr.set_ydata(tr.get_ydata()[:-1])
                        else:
                            continue

                tr.ydata = tr.ydata.astype(num.float)
                tr.ydata -= tr.ydata.mean()
                if self.tinc is not None:
                    win = num.hanning(tr.data_len())
                else:
                    win = num.ones(tr.data_len())

                tr.ydata *= win
                f, a = tr.spectrum(pad_to_pow2=True)
                a = num.abs(a)**2
                a *= tr.deltat * 2. / num.sum(win**2)
                a[0] /= 2.
                a[a.size//2] /= 2.

                if nslc not in by_nslc:
                    by_nslc[nslc] = []

                by_nslc[nslc].append((f, a))

        if not by_nslc:
            self.fail('No complete data windows could be exctracted for ' +
                      'given selection')

        fframe = self.figure_frame()
        fig = fframe.gcf()

        if self.join_stations:
            grouping = lambda k: (k[3],)
            labeling = lambda k: ' '.join(x for x in k[:-1] if x)
        else:
            grouping = lambda k: k
            labeling = lambda k: None

        group_keys = sorted(set(grouping(k) for k in by_nslc.keys()))

        p = None

        ncols = len(group_keys) // 5 + 1
        nrows = (len(group_keys)-1) // ncols + 1

        axes = []
        for i, group_key in enumerate(group_keys):
            p = fig.add_subplot(nrows, ncols, i+1, sharex=p, sharey=p)
            axes.append(p)

            legend = False
            for j, k in enumerate(sorted(by_nslc.keys())):
                color = to01(colors[j % len(colors)])
                color_trans1 = color + (0.5,)
                color_trans2 = color + (0.25,)

                group = by_nslc[k]
                if grouping(k) == group_key:
                    a_list = [a for (f, a) in group]
                    a = num.vstack(a_list)

                    if self.percentiles:
                        p10 = num.percentile(a, 10., axis=0)
                        p90 = num.percentile(a, 90., axis=0)
                        p.fill_between(
                            f[1:], p10[1:], p90[1:], color=color_trans1)

                    if self.minmax:
                        p0 = num.percentile(a, 0., axis=0)
                        p100 = num.percentile(a, 100., axis=0)
                        p.fill_between(
                            f[1:], p0[1:], p100[1:], color=color_trans2)

                    lab = labeling(k)
                    if self.mean:
                        mean = num.mean(a, axis=0)
                        p.plot(f[1:], mean[1:], label=lab, color=color)
                        if lab:
                            legend = True
                        lab = None

                    if self.logmean:
                        logmean = num.exp(num.mean(num.log(a), axis=0))
                        p.plot(f[1:], logmean[1:], label=lab, color=color)
                        if lab:
                            legend = True
                        lab = None

                    if self.median:
                        p50 = num.median(a, axis=0)
                        p.plot(f[1:], p50[1:], label=lab, color=color)
                        if lab:
                            legend = True
                        lab = None

            fmin = min(f[1] for (f, a) in group)
            fmax = max(f[-1] for (f, a) in group)
            if self.tinc is not None:
                fmin = max(fmin, 1.0/self.tinc)

            p.set_title(
                ' '.join(group_key), ha='right', va='top', x=0.99, y=0.9)
            p.grid()

            p.set_xscale('log')
            p.set_yscale('log')
            if i/ncols == (len(group_keys)-1)/ncols:
                p.set_xlabel('Frequency [Hz]')

            if i % ncols == 0:
                p.set_ylabel('PSD')

            p.set_xlim(fmin, fmax)

            if legend:
                p.legend(loc='lower left', prop=dict(size=9))

        for i, p in enumerate(axes):

            if i/ncols != (len(group_keys)-1)/ncols:
                for t in p.get_xticklabels():
                    t.set_visible(False)

            if i % ncols != 0:
                for t in p.get_yticklabels():
                    t.set_visible(False)
            else:
                tls = p.get_yticklabels()
                if len(tls) > 8:
                    for t in tls[1::2]:
                        t.set_visible(False)

        try:
            fig.tight_layout()
        except AttributeError:
            pass

        if self.save:
            fig.savefig(self.output_filename(dir='psd.pdf'))

        fig.canvas.draw()


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [PlotPSD()]

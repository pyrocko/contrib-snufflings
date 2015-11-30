import math
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
import numpy as num
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

from pyrocko import plot, util


def to01(c):
    return c[0]/255., c[1]/255., c[2]/255.


def desat(c, a):
    cmean = (c[0] + c[1] + c[2]) / 3.
    return tuple(cc*a + cmean*(1.0-a) for cc in c)


name_to_taper = {
    'Hanning': num.hanning,
    'Hamming': num.hamming,
    'Blackman': num.blackman,
    'Bartlett': num.bartlett}

cmap_colors = [plot.tango_colors[x] for x in [
    'skyblue1', 'chameleon1', 'butter1', 'orange1', 'scarletred1', 'plum3']]

name_to_cmap = {
    'spectro': LinearSegmentedColormap.from_list(
        'spectro', [desat(to01(c), 0.8) for c in cmap_colors])}


def get_cmap(name):
    if name in name_to_cmap:
        return name_to_cmap[name]
    else:
        return cm.get_cmap(name)


class Spectrogram(Snuffling):

    '''
    <html>
    <body>
    <h1>Plot spectrogram</h1>
    <p>Plots a basic spectrogram.</p>
    </body>
    </html>
    '''

    def setup(self):
        '''Customization of the snuffling.'''

        self.set_name('Spectrogram')
        self.add_parameter(
            Param('Window length [s]:', 'twin', 100, 0.1, 10000.))

        self.add_parameter(
            Param('Overlap [%]:', 'overlap', 75., 0., 99.))

        self.add_parameter(
            Switch('Save figure', 'save', False))

        self.add_parameter(
            Choice('Taper function', 'taper_name', 'Hanning',
                   ['Hanning', 'Hamming', 'Blackman', 'Bartlett']))

        self.add_parameter(
            Choice('Color scale', 'color_scale', 'log',
                   ['log', 'sqrt', 'lin']))

        self.add_parameter(
            Choice('Color table', 'ctb_name', 'spectro',
                   ['spectro', 'rainbow']))

        self.set_live_update(False)
        self._tapers = {}

    def get_taper(self, name, n):

        taper_key = (name, n)

        if taper_key not in self._tapers:
            self._tapers[taper_key] = name_to_taper[name](n)

        return self._tapers[taper_key]

    def call(self):
        '''Main work routine of the snuffling.'''

        by_nslc = {}
        tpad = self.twin * self.overlap/100. * 0.5
        tinc = self.twin - 2 * tpad
        times = []
        for traces in self.chopper_selected_traces(
                tinc=tinc, tpad=tpad, want_incomplete=False, fallback=True):

            for tr in traces:
                nslc = tr.nslc_id
                nwant = int(math.floor((tinc + 2*tpad) / tr.deltat))

                if nwant != tr.data_len():
                    if tr.data_len() == nwant + 1:
                        tr.set_ydata(tr.get_ydata()[:-1])
                    else:
                        continue

                tr.ydata = tr.ydata.astype(num.float)
                tr.ydata -= tr.ydata.mean()

                win = self.get_taper(self.taper_name, tr.data_len())
                tr.ydata *= win

                f, a = tr.spectrum(pad_to_pow2=True)
                df = f[1] - f[0]
                a = num.abs(a)**2
                a *= tr.deltat * 2. / (df*num.sum(win**2))
                a[0] /= 2.
                a[a.size/2] /= 2.

                if nslc not in by_nslc:
                    by_nslc[nslc] = []

                tmid = 0.5*(tr.tmax + tr.tmin)
                by_nslc[nslc].append((tmid, f, a))
                times.append(tmid)

        if not by_nslc:
            self.fail('No complete data windows could be exctracted for '
                      'given selection')

        fframe = self.figure_frame()
        fig = fframe.gcf()

        nslcs = sorted(by_nslc.keys())

        p = None

        ncols = len(nslcs) / 5 + 1
        nrows = (len(nslcs)-1) / ncols + 1

        tmin = min(times)
        tmax = max(times)
        nt = int(round((tmax - tmin) / tinc)) + 1
        t = num.linspace(tmin, tmax, nt)

        if (tmax - tmin) < 60:
            tref = util.day_start(tmin)
            tref += math.floor((tmin-tref) / 60.) * 60.
            t -= tref
            tunit = 's'
        elif (tmax - tmin) < 3600:
            tref = util.day_start(tmin)
            tref += math.floor((tmin-tref) / 3600.) * 3600.
            t -= tref
            t /= 60.
            tunit = 'min'
        else:
            tref = util.day_start(tmin)
            t -= tref
            t /= 3600.
            tunit = 'h'

        axes = []
        for i, nslc in enumerate(nslcs):
            p = fig.add_subplot(nrows, ncols, i+1, sharex=p, sharey=p)
            axes.append(p)
            group = by_nslc[nslc]
            f = group[0][1]

            nf = f.size
            a = num.zeros((nf, nt), dtype=num.float)
            a.fill(num.nan)
            for (t1, _, a1) in group:
                it = int(round((t1 - tmin) / tinc))
                if it < 0 or nt <= it:
                    continue

                a[:, it] = a1

            if self.color_scale == 'log':
                a = num.log(a)
                label = 'log PSD'
            elif self.color_scale == 'sqrt':
                a = num.sqrt(a)
                label = 'sqrt PSD'
            else:
                label = 'PSD'

            a = num.ma.masked_invalid(a)

            min_a = num.min(a)
            max_a = num.max(a)
            mean_a = num.mean(a)
            std_a = num.std(a)

            zmin = max(min_a, mean_a - 3.0 * std_a)
            zmax = min(max_a, mean_a + 3.0 * std_a)

            pcm = p.pcolormesh(t, f, a, cmap=get_cmap(self.ctb_name),
                               vmin=zmin, vmax=zmax)

            fmin = 2.0 / self.twin
            fmax = f[-1]

            p.set_title(
                '.'.join(x for x in nslc if x),
                ha='right',
                va='top',
                x=0.99,
                y=0.9)

            p.grid()

            p.set_yscale('log')

            divider = make_axes_locatable(p)
            cax = divider.append_axes('right', size='2%', pad=0.2)

            cbar = fig.colorbar(pcm, cax=cax)
            cbar.set_label(label)

            if i/ncols == (len(nslcs)-1)/ncols:
                p.set_xlabel('Time since %s [%s]' %
                             (util.time_to_str(tref, format='%Y-%m-%d %H:%M'),
                              tunit))

            if i % ncols == 0:
                p.set_ylabel('Frequency [Hz]')

            p.set_xlim(t[0], t[-1])
            p.set_ylim(fmin, fmax)

        for i, p in enumerate(axes):

            if i/ncols != (len(nslcs)-1)/ncols:
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
    return [Spectrogram()]

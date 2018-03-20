import logging
import numpy as num
import matplotlib.pyplot as plt
import copy

from matplotlib.collections import LineCollection
from pyrocko import util
from pyrocko.gui.snuffling import Snuffling, Choice, Switch


logger = logging.getLogger('pyrocko.gui.snuffling.particle_motion')


class ParticleMotion(Snuffling):
    '''
    Visualize Particle Motions

    Plot combinations of channels ending with one of: 'NEZ012' for traces
    selected either with an extended marker or which are visible in snuffler.
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        self.set_name('Particle Motion')
        self.add_parameter(
            Choice(
                'Colormap', 'cmap', 'viridis', ['viridis', 'coolwarm', 'cool',
                                                'bone', 'seismic'])
        )
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        viewer = self.get_viewer()
        nslc_ids = self.get_viewer().pile.nslc_ids

        nsl_ids = set([nslc_id[:3] for nslc_id in nslc_ids])
        figs = []
        for nsl_id in nsl_ids:
            x = []
            y = []
            z = []

            zlabel = ''
            xlabel = ''
            ylabel = ''

            def selector(tr):
                return util.match_nslc("%s.%s.%s.*" % nsl_id, tr.nslc_id)

            for traces in self.chopper_selected_traces(
                    fallback=True, trace_selector=selector):
                for tr in traces:
                    tr = tr.copy(data=True)
                    if viewer.highpass:
                        tr.highpass(4, viewer.highpass)
                    if viewer.lowpass:
                        tr.lowpass(4, viewer.lowpass)
                    time = tr.get_xdata()
                    net, sta, loc, cha = tr.nslc_id
                    if 'E' in cha or '1' in cha:
                        x = tr.ydata
                        xlabel = cha
                    elif 'N' in cha or '2' in cha:
                        y = tr.ydata
                        ylabel = cha
                    elif 'Z' in cha or '0' in cha:
                        z = tr.ydata
                        zlabel = cha
                    else:
                        self.warn(
                            'Did not find any character of "ENZ012" in %s' %
                            cha)
            is_zero = 0
            for val in [x, y, z]:
                is_zero += (len(val) == 0)
            if is_zero > 1:
                continue

            fig = self.pylab(get='figure')
            time -= min(time)
            axs = [None, None, None]
            for iax, (xi, yi, xilabel, yilabel) in enumerate([
                    (x, y, xlabel, ylabel),
                    (z, y, zlabel, ylabel),
                    (x, z, xlabel, zlabel),
            ]):
                ax = fig.add_subplot(
                    221 + iax, aspect='equal', sharex=axs[0], sharey=axs[0])

                if len(xi) == len(yi):
                    try:
                        mapable = ax.scatter(xi, yi, c=time, cmap=self.cmap)
                    except ValueError as e:
                        logger.error('nsl id: %s | %s' % (nsl_id, e))
                    ax.plot(xi, yi, c='grey', alpha=0.4)
                    ax.set_xlim(num.min(xi), num.max(xi))
                    ax.set_ylim(num.min(yi), num.max(yi))
                    ax.set_aspect('equal')
                ax.set_xlabel(xilabel)
                ax.set_ylabel(yilabel)
                ax.grid(True)
                if iax == 0:
                    ax.xaxis.set_label_position('top')
                    ax.xaxis.tick_top()
                elif iax == 1:
                    ax.xaxis.set_label_position('top')
                    ax.xaxis.tick_top()
                    ax.yaxis.set_label_position('right')
                    ax.yaxis.tick_right()

                axs[iax] = ax

                ax = fig.add_subplot(224)
                min_t = num.min(time)
                max_t = num.max(time)
                for i, (v, label) in enumerate(zip([x, y, z], [xlabel, ylabel, zlabel])):
                    v = copy.deepcopy(v)
                    mean_v = num.mean(v)
                    v -= mean_v
                    v /= num.max(v)
                    v += i
                    points = num.array([time, v]).T.reshape(-1, 1, 2)
                    segments = num.concatenate([points[:-1], points[1:]], axis=1)
                    lc = LineCollection(segments, cmap=plt.get_cmap(self.cmap))
                    lc.set_array(time)
                    ax.add_collection(lc)
                    ax.text(max_t, v[-1], label)
                    ax.set_xlabel('time [s]')

                ax.set_ylim(-1., 3.)
                ax.set_xlim(min_t, max_t)
                ax.spines['left'].set_color('none')
                ax.spines['top'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.yaxis.set_ticks([])

            fig.suptitle('.'.join(nsl_id))
            fig.subplots_adjust(hspace=0.05, wspace=0.05, bottom=0.07)
            figs.append(fig)

        for fig in figs:
            fig.canvas.draw()


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ParticleMotion()]

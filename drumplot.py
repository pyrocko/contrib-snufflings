import logging
import numpy as num
import matplotlib.pyplot as plt
import copy

from pyrocko import util, trace
from pyrocko.snuffling import Snuffling, Choice, Switch, Param

logger = logging.getLogger('pyrocko.snuffling.drumplot')


class DrumPlot(Snuffling):
    '''
    Drum Plot
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        self.set_name('Drum Plot')
        self.add_parameter(
            Param(
                'Gain', 'yscale', 2., 0.1, 100.)
        )
        self.add_parameter(
            Choice(
                'Pre-scale', 'prescale', 'max', ['max', 'std'])
        )
        self.add_parameter(
            Choice(
                'N minutes', 'xminutes', '15', ['1', '15', '30', '60'])
        )
        # self.add_parameter(
        #     Choice(
        #         'N hours', 'nhours', '24', ['24', '1'])
        # )
        self.add_parameter(Switch('Global common scale', 'scale_global', True))

        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        viewer = self.get_viewer()

        figs = {}
        fig_width_inch = viewer.width()
        npixel_hori = float(fig_width_inch*50)
        xminutes = int(self.xminutes)
        xseconds = xminutes * 60

        self.nhours = 24
        nrows = int(self.nhours) * 60 / xminutes
        ynormalizations = {}
        lines_data = {}

        for traces in self.chopper_selected_traces(tinc=60*60, fallback=True):
            for tr in traces:
                t0 = util.day_start(tr.tmin)
                key = (tr.nslc_id, t0)
                if key not in figs:
                    fig = self.pylab(get='figure')
                    ax = fig.add_subplot(111)
                    figs[key] = (fig, ax)
                    ynormalizations[key] = 0
                    lines_data[key] = []

                tr = tr.copy(data=True)
                ndecimate = int((xseconds/tr.deltat) / npixel_hori)
                tr.downsample(ndecimate)
                if self.prescale == 'max':
                    ynormalizations[key] = max(num.max(tr.ydata), ynormalizations[key])
                else:
                    ynormalizations[key] = max(num.std(tr.ydata), ynormalizations[key])

                if viewer.highpass:
                    tr.highpass(4, viewer.highpass)
                if viewer.lowpass and 1./tr.deltat>2.*viewer.lowpass:
                    tr.lowpass(4, viewer.lowpass)

                t = tr.get_xdata() - t0
                y = num.asarray(tr.get_ydata(), dtype=num.float)
                nskip = t / 3600.
                x = t % xseconds
                xdiff = num.diff(x) 
                itmp = num.where(num.logical_or(xdiff < 0, num.abs(xdiff-tr.deltat) > 1E-4))[0]
                indices = num.zeros(len(itmp)+2, dtype=num.int)
                indices[1:-1] = itmp
                indices[-1] = len(y)-1
                for i in range(len(indices)-1):
                    istart = indices[i] + 1
                    istop = indices[i+1]
                    lines_data[key].append(
                        (t0, x[istart: istop], y[istart: istop],
                         nskip[istart: istop])
                    )

        ynorm = None
        if self.scale_global:
            ynorm = num.max(ynormalizations.values())

        for key, lines in lines_data.items():
            for (t0, x, y, shifts) in lines:
                fig, ax = figs[key]
                ax.plot(
                    x/60.,
                    y/((ynorm or ynormalizations[key])/self.yscale) + shifts,
                    color='black')

                ax.set_title(util.tts(t0, format='%Y-%m-%d'))
        
        yticks = range(0, self.nhours+2, 2)
        xticks = range(0, xminutes+1, 1)
        for key, (fig, ax) in figs.items():
            ax.set_xlim(0, xminutes)
            ax.set_ylabel('Hour')
            ax.set_xlabel('Minute')
            ax.yaxis.set_ticks(yticks)
            ax.xaxis.set_ticks(xticks)
            ax.set_ylim(-0.1, 24.1)
            fig.canvas.draw()

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [DrumPlot()]

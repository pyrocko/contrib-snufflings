from pyrocko.gui.snuffling import Snuffling, Param, Switch
from collections import defaultdict
import copy
import numpy as num


class MarkerStack(Snuffling):

    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h2>Stack channels using selected markers</h2>
    Data segments are chopped using selected markers and padded as defined
    by 'Tmin padding' and 'Tmax padding'.
    </body>
    </html>
    '''
    def setup(self):
        self.set_name('Stack by markers')
        self.add_parameter(Param(
            'Tmin padding', 'tmin_pad', 1., 0., 100.))
        self.add_parameter(Param(
            'Tmax padding', 'tmax_pad', 1., 0., 100.))
        self.add_parameter(Switch('Normalize', 'normalize', False))
        self.add_parameter(Switch('Debug', 'debug', False))

    def call(self):
        self.cleanup()
        stencils = copy.deepcopy(self.get_selected_markers())
        for m in stencils:
            m.tmin -= self.tmin_pad
            m.tmax += self.tmax_pad

        t_first = min([m.tmin for m in stencils])

        pile = self.get_pile()
        stacks = defaultdict()
        nstacks = defaultdict()
        for m in stencils:
            for traces in pile.chopper(
                    tmin=m.tmin,
                    tmax=m.tmax,
                    trace_selector=lambda tr: m.match_nslc(tr.nslc_id)):

                for tr in traces:
                    nstacks.setdefault(tr.nslc_id, 0)
                    nstacks[tr.nslc_id] += 1

                    stack = stacks.get(tr.nslc_id, tr.copy(data=False))
                    stack.ydata = stack.ydata.astype(num.float)
                    stack.ydata -= num.mean(stack.ydata)
                    print('stacking %s ' % tr)
                    tr.shift(-(tr.tmin - t_first))

                    if self.debug:
                        print('adding trace:\n %s' % tr)
                        self.add_trace(tr.copy())
                        self.add_trace(tr)

                    stack.add(tr)
                    stacks[tr.nslc_id] = stack

        for nslc_id, s in stacks.items():
            if self.normalize:
                s.ydata /= nstacks[nslc_id]

            s.set_codes(location=s.location+'STACK')
            ymean = num.mean(s.ydata)
            s.ydata[0] = ymean
            s.ydata[-1] = ymean

        self.add_traces(stacks.values())


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''

    return [MarkerStack()]

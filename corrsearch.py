from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko.pile_viewer import Marker, EventMarker
from pyrocko import io, trace, model

class CorrsearchSnuffling(Snuffling):
    
    '''
    '''

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Cross Correlation Search')
        self.add_parameter(Param('Downsample to [Hz]', 'downsample', None, 0.1, 200., high_is_none=True))
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 1., 0.001, 50.))
        self.add_parameter(Param('Lowpass [Hz]', 'corner_lowpass', 1., 0.001, 50.))
        self.add_parameter(Switch('Apply to full dataset', 'apply_to_all', False))
        self.add_parameter(Choice('Normalization', 'normalization', 'Off', ('Off', 'Normal', 'Gliding')))
        self.add_parameter(Param('Threshold', 'threshold', 0.5, 0., 1.))
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        
        self.cleanup()
        
        period_highpass = 1./self.corner_highpass
        tpad = period_highpass
        
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
        for tr in pile.all(tmin=master_tmin, tmax=master_tmax, tpad=tpad):
            if self.downsample is not None:
                tr.downsample_to(1./self.downsample)
            tr.highpass(4, self.corner_highpass)
            tr.lowpass(4, self.corner_lowpass)
            tr.chop(tr.wmin, tr.wmax)
            masters[tr.nslc_id] = tr

        if self.apply_to_all:
            tmin, tmax = pile.get_tmin()+tpad, pile.get_tmax()
        else:
            tmin, tmax = self.get_viewer().get_time_range()
  
        tmaster = master_tmax-master_tmin
        tinc = min(20*tmaster, max(tmaster, tmax-tmin))

        for traces in pile.chopper(tmin=tmin, tmax=tmax, tinc=tinc, tpad=tmaster+tpad, want_incomplete=False):
            scc = None
            sccn = 0
            for b in traces:
                nslc = b.nslc_id
                if nslc in masters:
                    a = masters[nslc]
                    if self.downsample is not None:
                        b.downsample_to(1./self.downsample)
                    b.highpass(4, self.corner_highpass)
                    tr.lowpass(4, self.corner_lowpass)
                    normalization = {'Off': None, 'Normal': 'normal', 'Gliding': 'gliding'}[self.normalization]
                    c = trace.correlate(a,b, mode='valid', normalization=normalization)
                    c.shift(-c.tmin + b.tmin)
                    c.meta = { 'tabu' : True }
                    
                    if scc is None:
                        scc = c.copy()
                        scc.wmin = b.wmin
                        scc.wmax = b.wmax
                        scc.set_codes(network='', station='Sum Cross Correlation', location='', channel='')
                        sccn = 1

                    else:
                        scc.add(c)
                        sccn += 1
            
            if scc is not None:
                scc.ydata /= sccn
                scc.chop(scc.wmin, scc.wmax)

                markers = []
                for t, a in zip(*scc.peaks(self.threshold, tsearch=2./self.corner_highpass)):
                    m = EventMarker(model.Event(time=t, lat=0., lon=0., name='Event(%.2g)' % a))
                    markers.append(m)
                    
                if viewer:
                    self.add_traces([scc])
                    self.add_markers(markers)
                else:
                    io.save([scc], self.out_path, format='from_extension')
                 
                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ CorrsearchSnuffling() ]

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


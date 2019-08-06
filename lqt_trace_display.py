import logging
import numpy as np
import matplotlib.pyplot as plt
import copy

from matplotlib.collections import LineCollection
from pyrocko import util, trace
from pyrocko.gui.snuffling import Snuffling, Param, Choice, Switch


logger = logging.getLogger('pyrocko.gui.snuffling.lqt_trace_display')


class LQTRotation(Snuffling):
    '''
    Rotate 3-C N,E,Z-oriented waveforms into L-Q-T components, only for stations
    with *selected* markers in the Snuffler GUI.
    LQT rotation is based on the analysis of seismic signals within a time-window
    specified either (i) using an extended marker, or (i) using afixed time 
    interval in s.

    Plot combinations of channels ending with one of: 'NEZ012' for traces
    selected either with an extended marker or which are visible in snuffler.
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        self.set_name('LQT trace display')
        self.add_parameter(
            Param(
                'P-wave window length (s.)', 'winlen', 0.1, 0.01, 5.0, False, False, 
                False)
        )
        self.add_parameter(
            Choice(
                'Colormap', 'cmap', 'viridis', ['viridis', 'coolwarm', 'cool',
                                                'bone', 'seismic'])
        )
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        self.enable_pile_changed_notifications()
        viewer = self.get_viewer()  # Returns a PileViewerMain() instance
        # viewer.add_traces()
        winlen_s = self.get_settings()['winlen'] # P-wave window length in s.
        m = self.get_selected_markers()
        if len(m)==0:
            self.warn('No markers selected (click on markers, on use "a" key to select all.')
            pass


        for ma in m:
            nsl_ids = [ nlsc[:3] for nlsc in ma.nslc_ids ]  # (NET, LOC, STA, CHA) tuple
            nsl_id = nsl_ids[0]
            print(nsl_id)
            net, sta, loc = nsl_id

            x = []
            y = []
            z = []

            zlabel = ''
            xlabel = ''
            ylabel = ''

            def selector(tr):
                return util.match_nslc("%s.%s.%s.*" % nsl_id, tr.nslc_id)

            traces = viewer.pile.all(trace_selector=selector)
            for tr in traces:
                tr = tr.copy(data=True)
                dt = tr.deltat  # Sampling interval in s.
                #ns = np.ceil(winlen_s/dt)  # P-wave time-window length in nb of samples
                if viewer.highpass:
                    tr.highpass(4, viewer.highpass)
                if viewer.lowpass:
                    tr.lowpass(4, viewer.lowpass)
                time = tr.get_xdata()
                i_win = np.logical_and(time>=ma.tmin, time<=ma.tmin+winlen_s)
                net, sta, loc, cha = tr.nslc_id
                if (cha[2]=='X') or (cha[2]=='E') or (cha[2]=='1'):
                    x = np.reshape(tr.ydata, (1,-1))  # Restrict loading of traces to the P-wave time window
                    tr_x = tr
                    xlabel = cha
                elif (cha[2]=='Y') or (cha[2]=='N') or (cha[2]=='2'):
                    y = np.reshape(tr.ydata, (1,-1))  # idem
                    tr_y = tr
                    ylabel = cha
                elif (cha[2]=='Z') or (cha[2]=='0'):
                    z = np.reshape(tr.ydata, (1,-1)) # idem
                    tr_z = tr
                    zlabel = cha
                else:
                    self.warn(
                        'Did not find any character of "ENZ012" in %s' %
                        cha)
            is_zero = 0
            for val in [x, y, z]:
                is_zero += (len(val) == 0)
            if is_zero > 1:
                #self.warn('At least one component has no sample for station %s' % sta)
                continue

            xyz = np.concatenate((x,y,z), axis=0)
            xyz_win = xyz[:,i_win]
            # Compute the covariance matrix of the 3 components:
            cov = np.cov(xyz_win)
            print(cov)
            # Eigenvalue decomposition:
            vals, vecs = np.linalg.eig(cov)
            ip = np.argmax(vals)
            ldir = vecs[:,ip]  # P-wave direction (L)
            tdir = np.cross(ldir, [0,0,1])  # SH-wave direction (T)
            qdir = np.cross(ldir, tdir)  # SV-wave direction (Q): TODO Check the polarity !!
            # Rotate components into the LQT coordinate frame:
            l = np.dot(xyz.T, ldir)
            q = np.dot(xyz.T, qdir)
            t = np.dot(xyz.T, tdir)

            # Display rotated traces in a new panel:
            if False:
                ax = self.pylab(get='axes')
                lmax = np.max(l)
                qmax = np.max(l)
                tmax = np.max(l)
                maxi = max([lmax, qmax, tmax])
                ax.plot(time, l/maxi+3, label='L')
                ax.plot(time, q/maxi+2, label="Q")
                ax.plot(time, t/maxi+1, label="T")
                ax.set_xlabel('time [s]')
                ax.legend()

            # Add rotated traces to viewer: (can be later removed using self.cleanup() method)
            tr_l = trace.Trace(network=tr_x.network, station=tr_x.station,
                         location=tr_x.location, channel='__L',
                         tmin=tr_x.tmin, tmax=tr_x.tmax, 
                         deltat=tr_x.deltat, ydata=l)
            print(tr_l.name())
            tr_q = trace.Trace(network=tr_x.network, station=tr_x.station,
                         location=tr_x.location, channel='__Q',
                         tmin=tr_x.tmin, tmax=tr_x.tmax, 
                         deltat=tr_x.deltat, ydata=q)
            tr_t = trace.Trace(network=tr_x.network, station=tr_x.station,
                         location=tr_x.location, channel='__T',
                         tmin=tr_x.tmin, tmax=tr_x.tmax, 
                         deltat=tr_x.deltat, ydata=t)
            #self.add_traces([tr_l, tr_q, tr_t])
            trace.snuffle([tr_l, tr_q, tr_t], controls=False, markers=[ma])

    

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [LQTRotation()]

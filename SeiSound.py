from PyQt4.QtCore import QString, QThread, pyqtSlot, pyqtSignal, pyqtSlot, SIGNAL, QTimer
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
import pyrocko.trace as trace
from pyrocko.trace import CosFader
from scipy.io.wavfile import write
from pyrocko.gui_util import Marker
from scipy.signal import resample
import numpy as num
import os
import time


try:
    from PyQt4.phonon import Phonon
    no_phonon = False
except ImportError:
    no_phonon = True


class MarkerThread(QThread):
    def __init__(self, *args, **kwargs):
        self.viewer = kwargs.pop('viewer', None)
        self.media = kwargs.pop('media', None)
        QThread.__init__(self)

        self.marker = None
        self.timer = QTimer(self)
        self.connect(self.timer, SIGNAL('timeout()'), self.check_and_update)

    def handle_states(self, state):
        if state in [Phonon.PlayingState]:
            #self.goon = True
            self.check_and_update()

        if state == Phonon.LoadingState:
            if self.marker is None:
                self.marker = self.viewer.selected_markers()[0].copy()
                self.marker.tmax = self.marker.tmin
                self.tstart = self.marker.tmin
                self.viewer.add_marker(self.marker)

        if state == Phonon.StoppedState:
            #self.goon = False
            self.timer.stop()
            self.marker = None

    def check_and_update(self):
        tcurrent = self.media.currentTime()/1000.
        self.marker.tmin = self.tstart + tcurrent - self.tstart
        self.marker.tmax = self.tstart + tcurrent - self.tstart
        self.viewer.update()
        #if self.goon:
        self.timer.start(100)


class SeiSound(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>

    <h1 align="center">Export Traces to wav-file</h1>

    <body>
    <p>

    Mark a time range you would like to listen to with an extended marker and press 'Run'.<br>
    Use the scroll bar to *speed up* the recording by the chosen factor.<br>
    Click *Apply Main Control Filter* if you simply want to accept the main<br>
    control filter settings, or choose a different setting using the scroll bars.
    <p>
    Note, that this snuffling requires PyQt4's Phonon module to be installed.<br>

</body>
    '''
    def __init__(self, *args, **kwargs):
        Snuffling.__init__(self, *args, **kwargs)
        self.start = pyqtSignal()
        self.stop = pyqtSignal()

    def setup(self):
        self.set_name('Play/Save Audio')
        self.add_parameter(Param('speed up', 'speed_up', 1., 1., 10.))
        self.add_parameter(Choice('fps', 'fps', '9000',
                                  ('44100', '32000', '18000', '16000', '9000')))
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 0.1,
                                 0.001, 100., low_is_none=True))

        self.add_parameter(Param('Lowpass [Hz]', 'corner_lowpass', 10.,
                                 0.001, 100., high_is_none=True))

        self.add_parameter(Param('Fader [percentage]', 'tfade', 20, 0., 100.))
        self.add_parameter(Param('Volume', 'volume', 60., 0., 100.))
        self.add_trigger('Apply Main Control Filters', self.set_from_main)
        self.add_trigger('Export .wav file', self.export_wav)

        self.set_live_update(False)
        self._tmpdir = self.tempdir()
        self.marker = None

    def call(self):
        self.cleanup()
        if no_phonon:
            self.warn('Install pyqt4 phonon!\nCan only export wav files.')
            self.export_wav()
        else:
            self.play_phonon()

    def prepare_data(self):
        trange = self.get_selected_time_range()
        self.ttotal = float(trange[1]-trange[0])
        ntraces = 1
        nslc_ids = []
        for tr in self.chopper_selected_traces():
            if ntraces != 1:
                self.fail('Can only play one selected trace')
            t = tr[0].copy()
            nslc_ids.append(t.nslc_id)

            t.taper(CosFader(xfrac=self.tfade/100.))
            if self.corner_lowpass:
                t.lowpass(4, self.corner_lowpass)
            if self.corner_highpass:
                t.highpass(4, self.corner_highpass)
            data = t.get_ydata()
            print data
            ntraces += 1
        return nslc_ids, data

    def play_phonon(self):
        self.nslc_ids, data = self.prepare_data()
        tmpfile = os.path.join(self._tmpdir, 'phononquake.wav')
        self.export_wav(data=data, fn=tmpfile)
        output = Phonon.AudioOutput(parent=self._panel_parent)
        output.setVolume(self.volume/100)
        self.m_media = Phonon.MediaObject(self._panel_parent)
        Phonon.createPath(self.m_media, output)
        self.m_media.setCurrentSource(Phonon.MediaSource(tmpfile))
        self.viewer = self.get_viewer()
        self.marker_thread = MarkerThread(viewer=self.viewer, media=self.m_media)
        self.m_media.stateChanged.connect(self.marker_thread.handle_states)
        #self.m_media.connect(self.m_media, SIGNAL('stateChanged'), self.marker_thread.handle_states)
        #self.m_media.connect(self.m_media, SIGNAL('stateChanged'), self.handle_states)
        #self.m_media.stateChanged.connect(self.handle_states)
        self.m_media.play()

    def export_wav(self, data=None, fn=None):
        if fn is None:
            fn = self.output_filename()
        if data is None:
            nslc_ids, data = self.prepare_data()
        data = num.asarray(data, dtype=num.float)
        n = trace.nextpow2(len(data))
        data = num.append(data, num.zeros(n-len(data)))
        arg = int(int(self.fps)*self.ttotal/self.speed_up)
        data = resample(data, arg)
        scaled = num.int16(data/float(num.max(num.abs(data))) * 32767)
        write(fn, int(self.fps), scaled)

    def set_from_main(self):
        v = self.get_viewer()
        self.set_parameter('corner_highpass', v.highpass)
        self.set_parameter('corner_lowpass', v.lowpass)

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ SeiSound() ]


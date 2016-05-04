from PyQt4.QtCore import QThread, SIGNAL, QTimer
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
import pyrocko.trace as trace
from pyrocko.trace import CosFader
from scipy.io.wavfile import write
from scipy.signal import resample
import numpy as num
import os
import tempfile


try:
    from PyQt4.phonon import Phonon
    no_phonon = False
except ImportError:
    no_phonon = True


class MarkerThread(QThread):
    def __init__(self, *args, **kwargs):
        self.viewer = kwargs.pop('viewer', None)
        self.follow = kwargs.pop('follow', False)
        QThread.__init__(self)
        self.marker = None
        self.media = None
        self.speed_up = 1.
        self.timer = QTimer(self)
        self.timer.setInterval(100)
        self.connect(self.timer, SIGNAL('timeout()'), self.check_and_update)
        self.previous_state = Phonon.StoppedState
        self.t_stretch = 1.
        self.time_range = (0., 0.)

    def handle_states(self, state):
        if state == Phonon.PausedState:
            self.timer.stop()

        if state == Phonon.PlayingState:
            if self.timer.isActive() is False:
                self.timer.start()
            self.check_and_update()

        if state == Phonon.LoadingState:
            pass

        if state == Phonon.StoppedState:
            self.cleanup()

        self.previous_state = state

    def cleanup(self):
        self.viewer.remove_markers([self.marker])
        self.marker = None
        self.timer.stop()
        self.viewer.update()

    def check_and_update(self):
        if self.marker is None:
            self.marker = self.viewer.selected_markers()[0].copy()
            self.time_range = (self.marker.tmin, self.marker.tmax)
            self.viewer.add_marker(self.marker)
            self._factor = self.speed_up/(1-self.t_stretch)
            if self.speed_up<0.:
                self._start_at = self.time_range[1]
            else:
                self._start_at = self.time_range[0]

        tcurrent = self.media.currentTime()/1000.
        now = self._start_at + tcurrent * self._factor
        self.marker.tmin = now
        self.marker.tmax = now
        if self.follow:
            self.viewer.go_to_time(self.marker.tmin)
        self.viewer.update()


class SeiSound(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>

    <h1 align="center">Listen to your seismic records</h1>

    <body>
    <p>

    Mark a time range you would like to listen to with an extended marker and press 'Run'.<br>
    Use the scroll bar to <b>fast forward</b> the recording by the chosen factor.<br>
    Click <b>Apply Main Control Filter</b> if you want to adopt the main<br>
    control filter settings or choose a different setting using the scroll bars.
    <p>
    Note, that this snuffling requires PyQt4's Phonon module to be installed.<br>
    You can try installing it through your distributions package manager. On
    Ubuntu/Debian:
    <b><pre>
    sudo apt-get install python-qt4-phonon
    </pre></b>
    <p>
    Due to performance reasons of the resampling, the fast forward factor is
    rounded to one decimal place.
    </body>
    '''

    def setup(self):
        self.set_name('Play/Save Audio')
        self.add_parameter(Param('Fast Forward/Rewind', 'speed_up', 1., -20., 30.))
        self.add_parameter(Choice('fps', 'fps', '4000',
                                  ('44100', '32000', '18000', '16000', '9000',
                                   '4000')))
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 0.1,
                                 0.001, 100., low_is_none=True))

        self.add_parameter(Param('Lowpass [Hz]', 'corner_lowpass', 100.,
                                 0.001, 100., high_is_none=True))

        self.add_parameter(Param('Fader [percentage]', 'tfade', 5, 0., 50.))
        self.add_parameter(Param('Volume', 'volume', 60., 0., 100.))
        self.add_parameter(Switch('Follow', 'follow', False))
        self.add_trigger('Pause/Play', self.pause_play)
        self.add_trigger('Stop', self.stop_play)
        self.add_trigger('Apply Main Control Filters', self.set_from_main)
        self.add_trigger('Export .wav file', self.export_wav)

        self.set_live_update(False)
        self._tmpdir = self.tempdir()
        self.marker = None
        self.output = None
        self.marker_thread = None
        if not no_phonon:
            self.m_media = Phonon.MediaObject(self._panel_parent)
        self.no_phonon_warn = 'Install pyqt4 phonon!\nCan only export wav files.\nCheckout this snuffling\'s help.'

    def my_cleanup(self):
        if self.marker_thread is not None:
            self.marker_thread.cleanup()
        self.cleanup()

    def call(self):
        self.my_cleanup()
        self.viewer = self.get_viewer()
        if no_phonon:
            self.warn(self.no_phonon_warn)
            self.export_wav()
        else:
            self.marker_thread = MarkerThread(viewer=self.viewer, follow=self.follow)
            if self.m_media.state() == Phonon.PlayingState:
                self.m_media.stop()
                self.marker_thread.cleanup()
            self.play_phonon()

    def prepare_data(self):
        trange = self.get_selected_time_range()
        if num.abs(trange[0]-trange[1]) == 0:
            self.fail('no time range selected')

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
            ntraces += 1

        return nslc_ids, data

    def play_phonon(self):
        self.nslc_ids, data = self.prepare_data()
        tmpfile = tempfile.mkstemp(dir=self._tmpdir, suffix='.wav')[1]
        #tmpfile = '/tmp/test.wav'
        self.export_wav(data=data, fn=tmpfile)
        output = Phonon.AudioOutput(parent=self._panel_parent)
        output.setVolume(self.volume/100)
        Phonon.createPath(self.m_media, output)
        self.m_media.setCurrentSource(Phonon.MediaSource(tmpfile))
        self.marker_thread.media = self.m_media
        self.marker_thread.speed_up = speed_up=num.round(self.speed_up, 1) 
        self.m_media.stateChanged.connect(self.marker_thread.handle_states)
        self.m_media.play()

    def export_wav(self, data=None, fn=None):
        if fn is None:
            fn = self.output_filename()
        if data is None:
            nslc_ids, data = self.prepare_data()
        data = num.asarray(data, dtype=num.float)
        n = trace.nextpow2(len(data))-len(data)
        n_frac = float(n)/(n+len(data))
        data = num.append(data, num.zeros(n))
        fps = int(self.fps)
        arg = int(fps*self.ttotal/num.abs(num.round(self.speed_up, 1)))
        data = resample(data, arg)
        nnew = len(data)
        data = data[:-int(n_frac*nnew)]
        if self.speed_up<0.:
            data = data[::-1]
        data[0]  = 0.
        self.marker_thread.t_stretch = n_frac
        scaled = num.int16(data/float(num.max(num.abs(data))) * 32767)
        write(fn, fps, scaled)

    def set_from_main(self):
        v = self.get_viewer()
        self.set_parameter('corner_highpass', v.highpass)
        self.set_parameter('corner_lowpass', v.lowpass)

    def stop_play(self):
        if self.m_media is not None:
            self.m_media.stop()

    def pause_play(self):
        if no_phonon:
            self.fail(self.no_phonon_warn)
        if self.m_media is None:
            self.call()
            return 

        state = self.m_media.state()
        if state == Phonon.PlayingState:
            self.m_media.pause()
        elif state == Phonon.PausedState:
            self.m_media.play()
        elif state == Phonon.LoadingState:
            self.call()
        elif state == Phonon.StoppedState:
            self.call()
        else:
            print 'unexpected state. cleanup....'
            self.m_media.stop()

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ SeiSound() ]


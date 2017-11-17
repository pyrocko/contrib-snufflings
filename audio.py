from __future__ import print_function

from scipy.io.wavfile import write, read
from scipy.signal import resample
import numpy as num
import tempfile

from pyrocko.gui.qt_compat import qc
from pyrocko.snuffling import (Snuffling, Param, Choice, Switch,
                               NoTracesSelected)
import pyrocko.trace as trace
from pyrocko.trace import CosFader


try:
    from PyQt5.phonon import Phonon
    no_phonon = False
except ImportError:
    no_phonon = True


class MarkerThread(qc.QThread):
    def __init__(self, *args, **kwargs):
        self.viewer = kwargs.pop('viewer', None)
        self.follow = kwargs.pop('follow', False)
        qc.QThread.__init__(self)
        self.marker = None
        self.media = None
        self.speed_up = 1.
        self.timer = qc.QTimer(self)
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.check_and_update)
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
            if self.speed_up < 0.:
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

    Mark a time range you would like to listen to with an extended markerand press 'Run'.<br>
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
    </body>
    '''

    def setup(self):
        self.set_name('Play/Save Audio')
        self.add_parameter(Param('Fast Forward/Rewind', 'speed_up', 10., -20., 30.))
        self.add_parameter(Choice('fps', 'fps_choice', '16000',
                                  ('44100', '32000', '18000', '16000', '9000',
                                   '4000', 'keep original')))
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 0.001,
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
        self.add_trigger('Load .wav file', self.load_data)

        self.set_live_update(False)
        self._tmpdir = self.tempdir()
        self.output = None
        self.marker_thread = None
        if not no_phonon:
            self.m_media = Phonon.MediaObject(self._panel_parent)
        self.no_phonon_warn = 'Install pyqt4 phonon!\nCan only export wav files.\nCheckout this snuffling\'s help.'
        self.added_traces = []

    def my_cleanup(self):
        if self.marker_thread is not None:
            self.marker_thread.cleanup()
        self.cleanup()
        for tr in self.added_traces:
            self.add_trace(tr)

    def call(self):
        self.my_cleanup()
        self.viewer = self.get_viewer()
        try:
            trange = self.get_selected_time_range(fallback=False)
        except NoTracesSelected:
            self.fail('no time range selected')

        if no_phonon:
            self.warn(self.no_phonon_warn)
            self.export_wav()
        else:
            self.marker_thread = MarkerThread(viewer=self.viewer, follow=self.follow)
            if self.m_media.state() == Phonon.PlayingState:
                self.m_media.stop()
                self.marker_thread.cleanup()
            self.play_phonon()

    def load_data(self):
        fn = self.input_filename()
        sampling_rate, data = read(fn)
        for i, channel in enumerate(data.T):
            tr = trace.Trace(
                tmin=0., ydata=channel, deltat=1./sampling_rate, station='wav%s' % i)
            self.add_trace(tr)
            self.added_traces.append(tr)
        self.set_parameter('fps_choice', 'keep original')

    def prepare_data(self):
        trange = self.get_selected_time_range()
        self.ttotal = float(trange[1]-trange[0])
        ntraces = 1
        nslc_ids = []
        for tr in self.chopper_selected_traces():
            if ntraces != 1:
                self.fail('Can only play one selected trace')
            self.fps = 1./tr[0].deltat
            t = tr[0].copy()
            t.set_ydata(num.asarray(t.ydata-num.mean(t.ydata), dtype=num.float))
            nslc_ids.append(t.nslc_id)
            if self.corner_lowpass:
                t.lowpass(4, self.corner_lowpass)
            if self.corner_highpass:
                t.highpass(4, self.corner_highpass)
            t.taper(CosFader(xfrac=self.tfade/100.))
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
        self.marker_thread.speed_up = int(num.round(self.speed_up))
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
        if not self.fps_choice == 'keep original':
            fps = int(self.fps_choice)
            arg = int(fps*self.ttotal/num.abs(num.round(self.speed_up)))
            data = resample(data, arg)
            nnew = len(data)
            data = data[:-int(n_frac*nnew)]
            if self.speed_up<0.:
                data = data[::-1]
            data[0]  = 0.
            self.fps = fps

        scaled = num.int16(data/float(num.max(num.abs(data))) * 32767)
        write(fn, self.fps, scaled)

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
            print('unexpected state. cleanup....')
            self.m_media.stop()


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [ SeiSound() ]


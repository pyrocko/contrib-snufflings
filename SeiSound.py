from pyrocko.snuffling import Snuffling, Param, Choice, Switch
from scipy.io.wavfile import write
from scipy.signal import resample
import numpy as np

class SeiSound(Snuffling):
    '''
Export Traces to wav-file
=========================

Mark a time range you would like to listen to with an extended marker. Press 'Run' and open the exported file with a music player. 

The exported file named 'quake_sound.wav' will be stored in the directory from where you opened snuffler. 

Use the scroll bar to *speed up* the recording by the chosen factor.

A one second lasting wav file stored at 44100 Hz is about 5.1 MB large. Use the *fps dropdown menu* to reduce the resolution of your sound file, if desired. 
The lowpass and highpass filter of the main controls will be applied before saving the wav file.

Click *Apply Main Control Filter* if you simply want to accept the main control filter settings, or choose a different setting using the scroll bars. 

    '''
    
    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Export .wav')
        self.add_parameter(Param('speed up', 'speed_up', 1., 1., 10.))
        self.add_parameter(Choice('fps', 'fps', '44100', ('44100', '32000', '18000', '16000', '9000'))) 
        self.add_parameter(Param('Highpass [Hz]', 'corner_highpass', 1.0,
        0.001, 50., low_is_none=True))
        self.add_parameter(Param('Lowpass [Hz]', 'corner_lowpass', 4.0,
        0.001, 50., high_is_none=True))
        self.add_trigger('Apply Main Control Filter', self.set_from_main)

        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        self.cleanup()
        trange = self.get_selected_time_range()
        ttotal = trange[1]-trange[0]
        for tr in self.chopper_selected_traces():
            t = tr[0].copy()
            if self.corner_lowpass:
                t.lowpass(4, self.corner_lowpass)
            if self.corner_highpass:
                t.highpass(4, self.corner_highpass)
            data = t.get_ydata()
            data = resample(data, int(int(self.fps)*ttotal/float(self.speed_up)))
            scaled = np.int16(data/float(np.max(np.abs(data))) * 32767)

        write('quake_sound.wav', int(self.fps) , scaled)

    def set_from_main(self):
        v = self.get_viewer()
        self.set_parameter('corner_highpass', v.highpass)
        self.set_parameter('corner_lowpass', v.lowpass)

def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ SeiSound() ]


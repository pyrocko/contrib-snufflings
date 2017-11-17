from pyrocko import trace
import math
import numpy as num
from pyrocko.gui.snuffling import Param, Snuffling, Switch
import scipy.signal as S
import scipy.stats as SS


def detrend_data(x_axis,data):
    

    slope,offset,dummy,dummy2,dumy3 = SS.linregress(x_axis, data)
    
    detrended_data = S.detrend(data)
    
    return detrended_data, slope, offset


def retrend_data(x_axis,data,slope, offset):
    
    retrended_data = x_axis * slope + data + offset

    return retrended_data


class Notch(Snuffling):

    '''Apply gaussian notch filter to data in Snuffler
    
The filter is defined by a central frequency and a window width.  If possible,
harmonics are removed as well.  Filtering carried out in frequency domain.
Data trace is detrended, FFT calculated, respective notch coefficients
calculated and multiplied onto the FFT section. Then back-transform to time
domain, adding the underlying trendline.'''
    
    def setup(self):    
        '''Customization of Notch Snuffling.'''
        
        self.set_name('Notch Filter')

        self.add_parameter(Param('Center Frequency [Hz]', 
                'centerfreq',  50., 0.001, 1000.0, low_is_none=False))

        self.add_parameter(Param('Notch width (FWHM) [Hz]',  
                'notchwidth', 1., 0.001, 1000.0, low_is_none=False))

        self.add_parameter(Switch('Also filter harmonics',
                'filter_harmonics', False))

        self.set_have_post_process_hook(True)

        self.FFF = {}

    def call(self):
        # snuffling parameters have changed

        # release cached frequency responses
        self.FFF = {}

        # released cached processed traces in snuffler internals
        viewer = self.get_viewer()
        viewer.old_processed_traces = None

    def post_process_hook(self, traces):

        FFF = self.FFF

        for tr in traces:
            indi = num.arange(tr.data_len(), dtype=num.float)
            detr_data, m, b = detrend_data(indi, tr.get_ydata()) 

            tr.set_ydata(detr_data)

            freqs, FT_data = tr.spectrum(pad_to_pow2=True, tfade=None)

            fff_length = FT_data.size

            key = (tr.deltat, fff_length)

            if key not in FFF:
                FFF[key] = num.ones((fff_length), dtype=num.complex)
                fcenter  = self.centerfreq
                fwidth = self.notchwidth

                iharm = 1
                while iharm*fcenter - fwidth*2 < freqs[-1]:
                    FFF[key] *= GaussNotch(iharm*fcenter, fwidth).evaluate(freqs)
                    iharm += 1
                    if not self.filter_harmonics:
                        break

            filtered_data = num.fft.irfft(FT_data*FFF[key])[:tr.data_len()]
            
            retrended_data = retrend_data(indi, filtered_data, m, b)
            
            tr.set_ydata(retrended_data)

        return traces

class GaussNotch(trace.FrequencyResponse):
    def __init__(self, center, fwhm):
        self.center = center
        self.fwhm = fwhm

    def evaluate(self, freqs):
        denom = self.fwhm / (2.*math.sqrt(math.log(2.)))
        return 1.0 - num.exp(-((freqs-self.center)/denom)**2)

def __snufflings__():

    return [ Notch()]


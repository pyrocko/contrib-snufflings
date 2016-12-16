import numpy as num
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
import okada, numpy




class okadaforward(Snuffling):
    '''
    <html>
    <body>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Plot Okada forward modelling re-wrapped/h1>
    <p>
    Plot standard in C-Band
    </body>
    </html>
    '''
    def setup(self):
        self.set_name("Okada Forward")
        self.add_parameter(Param(
            'Dip', 't_dip', 45., 0., 90.))
        self.add_parameter(Param(
            'Strike', 't_strike', 130., 0.1, 360.))
        self.add_parameter(Param(
            'Rake', 't_rake', 120., -180., 180.))
        self.add_parameter(Param(
            'Slip', 't_slip', 1., 0.1, 20.))
        self.add_parameter(Param(
            'Depth of Top edge', 't_ztop', -1.5e3, -50e3, 0.))
        self.add_parameter(Param(
            'Depth of bottom edge', 't_zbot', -4e3, -50e3, 0.))
        self.add_parameter(Param(
            'Length', 't_length', 10e3, 1, 50e3))
        self.add_parameter(Param(
            'Grid extent', 't_ext', 25e3, 10e3, 200e3))
        self.add_parameter(Param(
            'X-shift', 't_xtrace', 0, 0., 200e3))
        self.add_parameter(Param(
            'Y-shfit', 't_ytrace', 0, 0., 200e3))

        self.add_trigger('L-Band', self.lband)
        self.add_trigger('Save Last Figure', self.save)
        self.set_live_update(False)
        self.fig = None

    def call(self):

        self.cleanup()
        viewer = self.get_viewer()
        

        los = 0.3815, 0.0843, 0.9205 # unit vector
        wavelength = .056 # meter C-Band
        extent = -self.t_ext, self.t_ext, -self.t_ext, self.t_ext # meter (xmin,xmax,ymin,ymax)
        
        fault = okada.OkadaSource(
          strike=self.t_strike, dip=self.t_dip, rake=self.t_strike, # degree
          slip=self.t_slip, # meter
          ztop=self.t_ztop, zbottom=self.t_zbot, length=self.t_length, # meter
          xtrace=self.t_xtrace, ytrace=self.t_ytrace ) # meter
        
        
        
        Y, X = numpy.meshgrid(
          numpy.linspace( extent[2], extent[3], 500 ),
          numpy.linspace( extent[0], extent[1], 500 ) )
                               
        XYZ = numpy.array([ X, Y, numpy.zeros_like(X) ]).T

        
        disp = fault.displacement( XYZ, poisson=.25 )
        

        
        disp_los = numpy.dot( disp, los )
        phase = ( numpy.mod( disp_los / ( .5 * wavelength ) * 2 + 1, 2 ) - 1 ) * numpy.pi
        
            
        from matplotlib import pylab as plt

        plt.ion()
        plt.subplot(111)
        plt.imshow( phase, extent=extent, cmap=plt.cm.jet, origin='lower' )
        plt.clim( [ -numpy.pi, numpy.pi ] )
        
        if not self._live_update:
            dx = numpy.array((-.5,.5)) * fault.length * numpy.sin( fault.strike * numpy.pi / 180 )
            dy = numpy.array((-.5,.5)) * fault.length * numpy.cos( fault.strike * numpy.pi / 180 )
            
            plt.plot( fault.xtrace + dx, fault.ytrace + dy, 'w-', linewidth=5, solid_capstyle='round' )
            plt.plot( fault.xtrace + dx, fault.ytrace + dy, 'k--', linewidth=2, dash_capstyle='round' )
            
        formatter = plt.FuncFormatter( lambda x, pos: '%dkm' % int( x / 1e3 ) if x else '0' )
        ax = plt.gca()
        ax.xaxis.set_major_formatter( formatter )
        ax.yaxis.set_major_formatter( formatter )
        
        plt.axis( extent )
        plt.grid()
        
        plt.show()

     
        

      #  self.fig.canvas.draw()
        if self._live_update:
            plt.canvas,show()
    def save(self):
        fn = self.output_filename('Select Filename', 'okada_forward_rewrapped.png')
        self.fig.savefig(fn, pad_inches=0.1, bbox_inches='tight', dpi=320)
        
    def lband(self):
        wavelength = .240 # meter L-Band



def __snufflings__():
    return [okadaforward()]
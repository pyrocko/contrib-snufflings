import numpy as num
from pyrocko.snuffling import Snuffling, Param, Choice, Switch
import okada, numpy
from pyrocko import io

from pyrocko import gf, util
from pyrocko.parimap import parimap
from pyrocko.gf import Range
from pyrocko import gf, moment_tensor as mtm, trace


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
    <h1 align="center">Plot rect. dislocation forward modelling /h1>
    <p>
    Plot a rectangular dislocation source in a elastic half-space.<br>
    Plots in a seperate window.    Standard figure is given in C-Band and the signal is re-wrapped with
    into LOS. <br>
        If <b>Auto-Run</b> is activated the figure is updated automatically when
    modifying a value on the panel. Note that the fault trace will not be drawn in this mode. <br>
    The save buttion will <b>Save</b> the displacment three traces (U,N,E)
    and will be saved to the current directory. <br>

    
    </body>
    </html>
    '''
    def setup(self):
        self.set_name("Geodetic Okada forward modelling")
        self.add_parameter(Param(
            'Dip [deg.]', 't_dip', 45., 0., 90.))
        self.add_parameter(Param(
            'Strike [deg.]', 't_strike', 130., 0.1, 360.))
        self.add_parameter(Param(
            'Rake [deg.]', 't_rake', 120., -180., 180.))
        self.add_parameter(Param(
            'Slip [m]', 't_slip', 1., 0.1, 20.))
        self.add_parameter(Param(
            'Depth of Top edge [m]', 't_ztop', -1.5e3, -50e3, 0.))
        self.add_parameter(Param(
            'Depth of bottom edge [m]', 't_zbot', -4e3, -50e3, 0.))
        self.add_parameter(Param(
            'Length [m]', 't_length', 10e3, 1, 50e3))
        self.add_parameter(Param(
            'Grid extent [m]', 't_ext', 25e3, 10e3, 200e3))
        self.add_parameter(Param(
            'X-shift of fault centre from 0 [m]', 't_xtrace', 0, 0., 200e3))
        self.add_parameter(Param(
            'Y-shfit of fault centre from 0 [m]', 't_ytrace', 0, 0., 200e3))
        self.add_parameter(Param(
            'Wavelength for rewrapping [m]', 't_wavelength', 0.056, 0., 0.325))
        self.add_parameter(Param(
            'LOS 1', 't_los1', 0.3815, 0., 1.))
        self.add_parameter(Param(
            'LOS 2', 't_los2', 0.0843, 0., 1.))
        self.add_parameter(Param(
            'LOS 3', 't_los3', 0.9205, 0., 1.))



        self.add_trigger('Save as displ. Traces', self.save)
        self.add_trigger('Save as LOS displ. Traces', self.savelos)
        self.set_live_update(False)
        self.fig = None

    def call(self):

        self.cleanup()
        viewer = self.get_viewer()
        

        los = self.t_los1, self.t_los2, self.t_los3 # unit vector
        wavelength = self.t_wavelength # meter C-Band
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
            plt.canvas.show()
    def save(self):
        fault = okada.OkadaSource(
          strike=self.t_strike, dip=self.t_dip, rake=self.t_strike, # degree
          slip=self.t_slip, # meter
          ztop=self.t_ztop, zbottom=self.t_zbot, length=self.t_length, # meter
          xtrace=self.t_xtrace, ytrace=self.t_ytrace ) # meter
        
        
        extent = -self.t_ext, self.t_ext, -self.t_ext, self.t_ext # meter (xmin,xmax,ymin,ymax)        
        Y, X = numpy.linspace( extent[2], extent[3] ),numpy.linspace( extent[0], extent[1]) 

                               
        XYZ = numpy.array([ X, Y, numpy.zeros_like(X) ]).T
        

        disp = fault.displacement( XYZ, poisson=.25 )
        tmint = util.str_to_time('1970-01-01 00:05:00.000')
        tr_U = trace.Trace(station='disp', channel='Z', deltat=0.5, tmin=tmint, ydata=disp[:,0])
        io.save([tr_U], 'up_displacement.mseed')
        tr_N = trace.Trace(station='disp', channel='N', deltat=0.5, tmin=tmint, ydata=disp[:,1])
        io.save([tr_N], 'north_displacement.mseed')
        tr_N = trace.Trace(station='disp', channel='E', deltat=0.5, tmin=tmint, ydata=disp[:,2])
        io.save([tr_N], 'east_displacement.mseed')
        
    def savelos(self):
        fault = okada.OkadaSource(
          strike=self.t_strike, dip=self.t_dip, rake=self.t_strike, # degree
          slip=self.t_slip, # meter
          ztop=self.t_ztop, zbottom=self.t_zbot, length=self.t_length, # meter
          xtrace=self.t_xtrace, ytrace=self.t_ytrace ) # meter
        
        
        extent = -self.t_ext, self.t_ext, -self.t_ext, self.t_ext # meter (xmin,xmax,ymin,ymax)        
        Y, X = numpy.linspace( extent[2], extent[3] ),numpy.linspace( extent[0], extent[1]) 

                               
        XYZ = numpy.array([ X, Y, numpy.zeros_like(X) ]).T
        los = self.t_los1, self.t_los2, self.t_los3 # unit vector
        

        disp = fault.displacement( XYZ, poisson=.25 )
        disp_los = numpy.dot( disp, los )
        tmint = util.str_to_time('1970-01-01 00:05:00.000')
        tr_U = trace.Trace(station='disp_los', channel='Z', deltat=0.5, tmin=tmint, ydata=disp_los[:,0])
        io.save([tr_U], 'up_displacement_los.mseed')
        tr_N = trace.Trace(station='disp_los', channel='N', deltat=0.5, tmin=tmint, ydata=disp_los[:,1])
        io.save([tr_N], 'north_displacement_los.mseed')
        tr_N = trace.Trace(station='disp_los', channel='E', deltat=0.5, tmin=tmint, ydata=disp_los[:,2])
        io.save([tr_N], 'east_displacement_los.mseed')
        
        

def __snufflings__():
    return [okadaforward()]
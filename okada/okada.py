import numpy, sys, os, ctypes, platform
from source import Source

libokada_file = 'libokada.' + { 'Linux': 'so'}[platform.system()]
libokada_path = os.path.join( os.path.dirname(__file__), libokada_file )
libokada = ctypes.cdll.LoadLibrary( libokada_path )


dsin = lambda x: numpy.sin( x * numpy.pi / 180. )
dcos = lambda x: numpy.cos( x * numpy.pi / 180. )
dtan = lambda x: numpy.tan( x * numpy.pi / 180. )


def tryfuncs( kwargs, *funcs ):
  for f in funcs:
    try:
      return f( **kwargs )
    except TypeError:
      pass
  raise Exception, 'not supported: ' + ', '.join( kwargs.keys() )


class OkadaSource( Source ):

  def __init__( self, **kwargs ):
    strike, dip, length, width, xbottom, ybottom, zbottom, geom_args = tryfuncs( kwargs,
      lambda strike, dip, length, width, xtop, ytop, ztop, **leftover: ( strike, dip, length, width,
        xtop + dcos(dip) * dcos(strike) * width, ytop - dcos(dip) * dsin(strike) * width, ztop - dsin(dip) * width, leftover ),
      lambda strike, dip, length, width, bottom, **leftover: ( strike, dip, length, width ) + tuple(bottom) + ( leftover, ),
      lambda strike, dip, length, width, xbottom, ybottom, zbottom, **leftover: ( strike, dip, length, width,
        xbottom, ybottom, zbottom, leftover ),
      lambda strike, dip, zbottom, ztop, length, xtrace, ytrace, **leftover: ( strike, dip, length, (ztop-zbottom) / dsin(dip),
        xtrace - zbottom * dcos(strike) / dtan(dip), ytrace + zbottom * dsin(strike) / dtan(dip), zbottom, leftover ),
    )
    strikeslip, dipslip, opening, leftover = tryfuncs( geom_args,
      lambda strikeslip, dipslip, opening=0., **leftover: ( strikeslip, dipslip, opening, leftover ),
      lambda slip, rake, opening=0., **leftover: ( slip * dcos(rake), slip * dsin(rake), opening, leftover ),
    )
    assert not leftover, 'leftover arguments: %s' % ', '.join( leftover.keys() )
    self.params = numpy.array( [ strike, dip, length, width, xbottom, ybottom, zbottom, strikeslip, dipslip, opening ], dtype=float )
    for key, value in kwargs.items():
      numpy.testing.assert_almost_equal( value, getattr(self,key) )

  strike     = property( lambda self: self.params[0] )
  dip        = property( lambda self: self.params[1] )
  length     = property( lambda self: self.params[2] )
  width      = property( lambda self: self.params[3] )
  xbottom    = property( lambda self: self.params[4] )
  ybottom    = property( lambda self: self.params[5] )
  zbottom    = property( lambda self: self.params[6] )
  strikeslip = property( lambda self: self.params[7] )
  dipslip    = property( lambda self: self.params[8] )
  opening    = property( lambda self: self.params[9] )
  bottom     = property( lambda self: self.params[4:7] )
  area       = property( lambda self: self.length * self.width )
  slip       = property( lambda self: numpy.hypot( self.strikeslip, self.dipslip ) )
  rake       = property( lambda self: numpy.arctan2( self.dipslip, self.strikeslip ) * 180 / numpy.pi )
  top        = property( lambda self: self.bottom + self.width * self.dipvec )
  xtop       = property( lambda self: self.top[0] )
  ytop       = property( lambda self: self.top[1] )
  ztop       = property( lambda self: self.top[2] )
  center     = property( lambda self: self.bottom + .5 * self.width * self.dipvec )
  xtrace     = property( lambda self: self.xbottom + self.zbottom * dcos(self.strike) / dtan(self.dip) )
  ytrace     = property( lambda self: self.ybottom - self.zbottom * dsin(self.strike) / dtan(self.dip) )
  dipvec     = property( lambda self: numpy.array( [ -dcos(self.dip) * dcos(self.strike), dcos(self.dip) * dsin(self.strike), dsin(self.dip) ] ) )
  strikevec  = property( lambda self: numpy.array( [ dsin(self.strike), dcos(self.strike), 0. ] ) )
  openvec    = property( lambda self: numpy.array( [ dcos(self.strike) * dsin(self.dip), -dsin(self.strike) * dsin(self.dip), dcos(self.dip) ] ) )
  slipvec    = property( lambda self: self.strikeslip * self.strikevec + self.dipslip * self.dipvec + self.opening * self.openvec )

  def patches( self, n, m ):
    length = self.length / float(n)
    width = self.width / float(m)
    patches = []
    for i in range(n):
      for j in range(m):
        bottom = self.bottom + self.strikevec*((i+.5-.5*n)*length) + self.dipvec*j*width
        patch = OkadaSource(
          strike=self.strike, dip=self.dip,
          length=length, width=width,
          bottom=bottom,
          strikeslip=self.strikeslip, dipslip=self.dipslip )
        patches.append( patch )
    return patches

  @property
  def corners( self ):
    return self.bottom + numpy.array(
      [[ -.5 * self.length * self.strikevec, -.5 * self.length * self.strikevec + self.width * self.dipvec ],
       [ +.5 * self.length * self.strikevec, +.5 * self.length * self.strikevec + self.width * self.dipvec ]] )

  def _call( self, func, xyz, poisson, *shape ):
    xyz = numpy.ascontiguousarray( xyz, dtype=float )
    assert xyz.shape[-1] == 3
    out = numpy.empty( xyz.shape + shape, dtype=numpy.double )
    func( out.ctypes, self.params.ctypes, xyz.ctypes, ctypes.c_double(poisson), ctypes.c_int(xyz.size//3) )
    return out
  
  def displacement( self, xyz, poisson ):
    # void get_displacements( double *out, OkadaSource *src, double *where, double poisson, int count );
    return self._call( libokada.get_displacements, xyz, poisson )
  
  def gradient( self, xyz, poisson ):
    # void get_gradients( double *out, OkadaSource *src, double *where, double poisson, int count );
    return self._call( libokada.get_gradients, xyz, poisson, 3 ).swapaxes(-1,-2)

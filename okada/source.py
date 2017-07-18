import numpy


class Source(object):

  def strain( self, xyz, poisson ):
    grad = self.gradient( xyz, poisson )
    return .5 * ( grad + grad.swapaxes(-1,-2) )

  def stress( self, xyz, poisson, young ):
    lmbda = (poisson*young)/float((1+poisson)*(1-2*poisson))
    mu = young/float(2*(1+poisson))
    strain = self.strain( xyz, poisson )
    stress = (2*mu) * strain
    diag(stress)[...] += lmbda * numpy.trace( strain, axis1=-2, axis2=-1 )
    return stress

  def __add__( self, other ):
    if other is 0: # to allow sum
      return self
    assert isinstance( other, Source )
    return AddSource( self, other )

  def __radd__( self, other ):
    return self.__add__( other )

  def __mul__( self, other ):
    assert isinstance( other, (int,float) )
    if other == 1:
      return self
    return ScaleSource( self, other )

  def __rmul__( self, other ):
    return self.__mul__( other )


class AddSource(Source):

  def __init__( self, source1, source2 ):
    assert isinstance( source1, Source )
    assert isinstance( source2, Source )
    self.source1 = source1
    self.source2 = source2

  def displacement( self, xyz, poisson ):
    return self.source1.displacement( xyz, poisson ) \
         + self.source2.displacement( xyz, poisson )

  def gradient( self, xyz, poisson ):
    return self.source1.gradient( xyz, poisson ) \
         + self.source2.gradient( xyz, poisson )


class ScaleSource(Source):

  def __init__( self, source, scale ):
    assert isinstance( source, Source )
    assert isinstance( scale, (int,float) )
    self.source = source
    self.scale = scale

  def displacement( self, xyz, poisson ):
    return self.scale * self.source.displacement( xyz, poisson )

  def gradient( self, xyz, poisson ):
    return self.scale * self.source.gradient( xyz, poisson )

  def __mul__( self, other ):
    return self.source.__mul__( self.scale * scale )


def diag(A):

  assert A.shape[-1] == A.shape[-2]
  return numpy.lib.stride_tricks.as_strided( A, shape=A.shape[:-1], strides=A.strides[:-2]+(A.strides[-2]+A.strides[-1],) )

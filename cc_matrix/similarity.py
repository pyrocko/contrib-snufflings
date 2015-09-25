from pyrocko import trace, model
from pyrocko.guts import *
from pyrocko.gf.seismosizer import SeismosizerTrace, Target

class Similarity(Object):
    '''CC result of one event pair at one target.'''
    ievent = Int.T(help='Event identifier')
    jevent = Int.T(help='Event identifier')
    itarget = Int.T(help='Target identifier')
    cross_correlation = Float.T()
    cross_correlation_trace = SeismosizerTrace.T(optional=True)
    relative_amplitude = Float.T(help='Relative amplitude with reference to ievent.')
    time_lag = Float.T(help='Time lag at maximum of cross correlation')

class SimilarityMatrix(Object):
    ''' A container class to store a list of :py:class:`Similarity` instances
    and how they have been calculated.'''
    events = List.T(model.Event.T())
    targets = List.T(Target.T())
    similarities = List.T(Similarity.T())
    filters = List.T(trace.FrequencyResponse.T())
    padding = Float.T()
    windowing_method = String.T()
    vmin = Float.T()
    vmax = Float.T()


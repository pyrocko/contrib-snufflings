from pyrocko import trace, model
from pyrocko.guts import *
from pyrocko.gf.seismosizer import SeismosizerTrace, Target

class Similarity(Object):
    ievent = Int.T()
    jevent = Int.T()
    itarget = Int.T()
    cross_correlation = Float.T()
    cross_correlation_trace = SeismosizerTrace.T()
    relative_amplitude = Float.T()

class SimilarityMatrix(Object):
    events = List.T(model.Event.T())
    targets = List.T(Target.T())
    similarities = List.T(Similarity.T())
    filters = List.T(trace.FrequencyResponse.T())
    padding = Float.T()



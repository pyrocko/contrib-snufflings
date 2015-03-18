from similarity import SimilarityMatrix, Similarity
from pyrocko.guts import *
from pyrocko import trace
import sys


if __name__=='__main__':
    fn = sys.argv[1]
    matrix = load(filename=fn)
    traces = []
    for sim in matrix.similarities:
        traces.append(sim.cross_correlation_trace.pyrocko_trace())

    trace.snuffle(traces, events=matrix.events) 


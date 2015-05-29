from pyrocko.snuffling import Snuffling, Param
import sys 
import vtk 
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as num 
from pyrocko import model
from pyrocko import orthodrome as ortho
import os
from vtk_focmec import *

class VtkTest(Snuffling):
    def setup(self):
        self.set_name('VTK')
        self.add_parameter(Param('dummy', 'dummy', 0, 0, 100))


    def call(self):
        self.cleanup()
        viewer = self.get_viewer()
        markers = self.get_selected_event_markers()
        events = [m.get_event() for m in markers]
        #tmin, tmax = self.get_selected_time_range(fallback=True)
        #event_markers = filter(lambda x: x.tmin>=tmin and x.tmax<=tmax,
        #                       viewer.markers)

        #event_markers = filter(lambda x: isinstance(x, EventMarker), event_markers)
        #if self.maxd:
        #    event_markers = filter(lambda x: distance(self, x._event)<=self.maxd*km,
        #                       event_markers)

        #if event_markers==[]:
        #    self.fail('No events in selected area found')


        self.frame = QtGui.QFrame(self._panel_parent)
 
        self.vl = QtGui.QGridLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vtkWidget.setMinimumSize(900, 880)
        print dir(self.vtkWidget)
        self.vl.addWidget(self.vtkWidget)
 
        self.ren = vtk.vtkRenderer()
        print dir(self.ren)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        #self.ren.AddActor(actor)
        actors = [] 
        normals, centers, colors = read_data(events=events)
        colors = [(0,1,0)]*len(colors)
        kwargs = {"centers": centers, 'normals':normals, "colors":colors, "return_pdm":True}
        actor1, apd = make_polydata_actor(**kwargs)
        actors.append(actor1)

        setup_renderer(self.ren, actors, bboxpolydata=apd)


        self.ren.ResetCamera()

        #self.ren.SetViewport(0.,0.,1.,1.)
        self._panel_parent.add_tab('pressure', self.frame)
        print self.ren.GetSize()

        self.iren.Initialize()


def __snufflings__():
    return [VtkTest()]

